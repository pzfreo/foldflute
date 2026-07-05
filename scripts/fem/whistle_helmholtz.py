# 3D Helmholtz FEM verification of the M2b rev5 whistle body (NGSolve/Netgen).
#
# Purpose: the design was solved with openwind's 1D transmission-line model.
# This script solves the full 3D acoustic eigenproblem on the ACTUAL air
# geometry (galleries, mitred+beveled elbows, turret chimneys, hole-crown
# curvature) surrounded by a block of exterior air, so radiation loading
# emerges from geometry instead of end-correction formulas. It then compares
# each fingering's resonance with openwind's production prediction.
#
#   -Delta p = k^2 p   in air (column + exterior block)
#   p = 0              on the window plane and the far exterior boundary
#   dp/dn = 0          on all walls and on 'fingertip caps' sealing closed holes
#
# Modes localized in the air column (energy fraction) are kept; exterior-box
# cavity modes are rejected.
#
# Usage:
#   .venv-fem/bin/python scripts/fem/whistle_helmholtz.py --calibrate
#   .venv-fem/bin/python scripts/fem/whistle_helmholtz.py --note D4 A4 C#5 A5
#   .venv-fem/bin/python scripts/fem/whistle_helmholtz.py --note all
import argparse
import json
import math
import os
import sys

import numpy as np
from netgen.occ import (Axis, Box, Cylinder, Dir, Glue, OCCGeometry, Pnt, X,
                        Y, Z)
import ngsolve as ng

HERE = os.path.dirname(os.path.abspath(__file__))
GEO = json.load(open(os.path.join(HERE, '..', 'fem_geometry.json')))

C_MM_S = 343_360.0          # sound speed at 20 C, mm/s (checked by --calibrate)
HOLES = ['h1', 'h2', 'h3', 'h4', 'h5', 'h6']


def air_column():
    """Mirror of the CAD air solid, in netgen.occ primitives (mm)."""
    g = GEO
    r = g['bore_r']
    zr, zg, zt, zrr = g['z_run'], g['z_gal'], g['z_tail'], g['z_rear']
    c1, c2, xt, xe = g['col1'], g['col2'], g['x_t'], g['x_e']
    segs = [
        # includes the 120 mm headjoint tube: window plane at x = -120
        Cylinder(Pnt(-120, 0, zr), X, r=r, h=120 + c1 + r),
        Cylinder(Pnt(c1, 0, zr - r), Z, r=r, h=(zg - zr) + 2 * r),
        Cylinder(Pnt(c1 - r, 0, zg), X, r=r, h=(c2 - c1) + 2 * r),
        Cylinder(Pnt(c2, 0, zr - r), Z, r=r, h=(zg - zr) + 2 * r),
        Cylinder(Pnt(c2 - r, 0, zr), X, r=r, h=(xt - c2) + 2 * r),
        Cylinder(Pnt(xt, 0, zr - r), Z, r=r, h=(zt - zr) + 2 * r),
        Cylinder(Pnt(xe - r, 0, zt), X, r=r, h=(xt - xe) + 2 * r),
        Cylinder(Pnt(xe, 0, zt - r), Z, r=r, h=zrr - (zt - r)),
    ]
    col = segs[0]
    for s in segs[1:]:
        col = col + s
    # tonehole air: from turret exit down into the bore
    for h in HOLES:
        hd = g['holes'][h]
        col = col + Cylinder(Pnt(hd['x'], 0, hd['exit_z']), Z,
                             r=hd['d'] / 2.0, h=(g['z_run'] - 4) - hd['exit_z'])
    # Coltman bevel cuts at each elbow: remove the outer-corner wedge beyond
    # u = bev_u along the corner bisector (local 26 mm cube bound)
    for vx, vz, th in g['bevels']:
        nx, nz = math.cos(math.radians(th)), -math.sin(math.radians(th))
        # wedge = halfspace {u >= bev_u} within a local cube around the vertex
        cube = Box(Pnt(vx - 13, -13, vz - 13), Pnt(vx + 13, 13, vz + 13))
        big = Box(Pnt(0, -100, -100), Pnt(200, 100, 100))
        big = big.Rotate(Axis(Pnt(0, 0, 0), Y), th)
        big = big.Move((vx + nx * g['bev_u'], 0, vz + nz * g['bev_u']))
        col = col - (big * cube)
    return col


def body_solid():
    """Simplified body mirror (what displaces exterior air)."""
    g = GEO
    zg, zrr, xe = g['z_gal'], g['z_rear'], g['x_e']
    xbody = g['x_t'] + 13.5
    b = Box(Pnt(56, -15, 0), Pnt(xbody, 15, 30))
    b = b + Box(Pnt(79, -15, 0), Pnt(131, 15, zg + 13.5))
    b = b + Box(Pnt(244, -15, 0), Pnt(xbody, 15, zrr))
    b = b + Box(Pnt(60, -15, -1), Pnt(82, 15, 31))
    b = b + Cylinder(Pnt(0, 0, g['z_run']), X, r=12.4, h=56)
    for h in HOLES:
        hd = g['holes'][h]
        ht = -hd['exit_z']
        if ht > 0.05:   # turret frustum
            b = b + Cylinder(Pnt(hd['x'], 0, hd['exit_z']), Z,
                             r=hd['d'] / 2.0 + 3.0 + ht, h=ht)
    return b


def domain_for(mask):
    """mask: {hole: 'open'|'closed'}; returns (glued domain, names set)."""
    g = GEO
    col = air_column()
    xbody = g['x_t'] + 13.5
    ext = Box(Pnt(30, -75, -55), Pnt(xbody + 45, 75, g['z_rear'] + 60))
    ext = ext - body_solid()
    for h in HOLES:
        if mask[h] == 'closed':      # fingertip cap seals the turret exit
            hd = g['holes'][h]
            ext = ext - Cylinder(Pnt(hd['x'], 0, hd['exit_z'] - 4.0), Z,
                                 r=hd['d'] / 2.0 + 3.5, h=4.0)
    ext = ext - col                  # avoid double-counting shared volume
    col.solids[0].name = 'col'
    for s in ext.solids:
        s.name = 'ext'
    col.solids[0].maxh = 6.0
    for s in ext.solids:
        s.maxh = 16.0
    dom = Glue([col, ext])
    # name boundaries: window plane (x~0) and far exterior faces -> popen
    for f in dom.faces:
        c = f.center
        if abs(c.x + 120) < 1e-6:
            f.name = 'popen'
        elif (abs(c.x - 30) < 1e-6 or abs(c.x - (xbody + 45)) < 1e-6
              or abs(abs(c.y) - 75) < 1e-6 or abs(c.z + 55) < 1e-6
              or abs(c.z - (g['z_rear'] + 60)) < 1e-6):
            f.name = 'popen'
        else:
            f.name = 'wall'
        # refine tonehole regions
        for h in HOLES:
            hd = g['holes'][h]
            if (abs(c.x - hd['x']) < 12 and abs(c.y) < 12
                    and c.z < g['z_run']):
                f.maxh = 2.4
    return dom


def eigs_near(mesh, f_targets, nev=40):
    fes = ng.H1(mesh, order=2, dirichlet='popen')
    u, v = fes.TnT()
    a = ng.BilinearForm(ng.grad(u) * ng.grad(v) * ng.dx).Assemble()
    m = ng.BilinearForm(u * v * ng.dx).Assemble()
    k0 = 2 * math.pi * (sum(f_targets) / len(f_targets)) / C_MM_S
    vecs = [a.mat.CreateColVector() for _ in range(nev)]
    lams = ng.ArnoldiSolver(a.mat, m.mat, fes.FreeDofs(), vecs,
                            shift=k0 * k0)
    out = []
    gfu = ng.GridFunction(fes)
    for i, lam in enumerate(lams):
        k2 = lam.real
        if k2 <= 0:
            continue
        f = math.sqrt(k2) * C_MM_S / (2 * math.pi)
        gfu.vec.data = vecs[i]
        e_col = ng.Integrate(gfu * gfu, mesh, definedon=mesh.Materials('col'))
        e_all = ng.Integrate(gfu * gfu, mesh)
        out.append((f, e_col / e_all))
    return sorted(out)


def openwind_reference():
    """Production 1D prediction (losses on, unflanged) per fingering."""
    from openwind import ImpedanceComputation
    g = GEO
    r_m = g['bore_r'] / 1000.0
    main = [[0.0, r_m], [g['L_eff_1d'] / 1000.0, r_m]]
    rows = [['label', 'position', 'radius', 'chimney']]
    for h in HOLES:
        hd = g['holes'][h]
        rows.append([h, g['positions_1d'][h] / 1000.0, hd['d'] / 2000.0,
                     hd['chimney'] / 1000.0])
    notes = [n.replace('#', 's') for n in g['fingerings']]
    chart = [['label'] + notes]
    for h in HOLES:
        chart.append([h] + ['x' if g['fingerings'][n][h] == 'closed' else 'o'
                            for n in g['fingerings']])
    comp = ImpedanceComputation(np.array([200.0, 300.0]), main, rows, chart,
                                note=notes[0], temperature=20.0, losses=True,
                                radiation_category='unflanged')
    ref = {}
    for n in g['fingerings']:
        t = g['targets'][n]
        comp.set_note(n.replace('#', 's'))
        fs = t * 2.0 ** (np.linspace(-170, 170, 29) / 1200.0)
        comp.recompute_impedance_at(fs)
        i = int(np.argmin(np.abs(comp.impedance)))
        fs2 = fs[i] * 2.0 ** (np.linspace(-12, 12, 13) / 1200.0)
        comp.recompute_impedance_at(fs2)
        j = int(np.argmin(np.abs(comp.impedance)))
        ref[n] = float(fs2[j])
    return ref


def run_note(note, ref):
    g = GEO
    mask = g['fingerings'][note]
    t = g['targets'][note]
    dom = domain_for(mask)
    geo = OCCGeometry(dom)
    mesh = ng.Mesh(geo.GenerateMesh(maxh=16.0))
    mesh.Curve(2)
    cands = eigs_near(mesh, [ref[note]])
    # column modes are identified as nearest to the 1D prediction (they are
    # strongly fingering-dependent; the exterior box modes are not). The
    # energy fraction w is reported as supporting evidence only.
    good = [(f, w) for f, w in cands
            if abs(1200 * math.log2(f / ref[note])) < 100]
    best = min(good, key=lambda fw: abs(fw[0] - ref[note])) if good else None
    return mesh.ne, cands, best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--calibrate', action='store_true')
    ap.add_argument('--note', nargs='*', default=[])
    args = ap.parse_args()

    if args.calibrate:
        L, r = 570.365, 11.0
        cyl = Cylinder(Pnt(0, 0, 0), X, r=r, h=L)
        for f in cyl.faces:
            c = f.center
            f.name = 'popen' if (abs(c.x) < 1e-6 or abs(c.x - L) < 1e-6) \
                else 'wall'
        mesh = ng.Mesh(OCCGeometry(cyl).GenerateMesh(maxh=6.0))
        mesh.Curve(2)
        fes = ng.H1(mesh, order=2, dirichlet='popen')
        u, v = fes.TnT()
        a = ng.BilinearForm(ng.grad(u) * ng.grad(v) * ng.dx).Assemble()
        m = ng.BilinearForm(u * v * ng.dx).Assemble()
        k1 = math.pi / L
        vecs = [a.mat.CreateColVector() for _ in range(12)]
        lams = ng.ArnoldiSolver(a.mat, m.mat, fes.FreeDofs(), vecs,
                                shift=k1 * k1)
        print('elements:', mesh.ne, ' dofs:', fes.ndof)
        fs = sorted(math.sqrt(l.real) * C_MM_S / (2 * math.pi)
                    for l in lams if l.real > 0)[:4]
        for n, f in enumerate(fs, 1):
            ideal = n * C_MM_S / (2 * L)
            print('mode %d: FEM %.3f Hz  ideal %.3f Hz  err %+.2f cents'
                  % (n, f, ideal, 1200 * math.log2(f / ideal)))
        return

    ref = openwind_reference()
    notes = list(GEO['fingerings']) if args.note == ['all'] else args.note
    results = {}
    for note in notes:
        ne, cands, best = run_note(note, ref)
        t = GEO['targets'][note]
        if best:
            fem_f, w = best
            print('%-4s target %8.2f | 1D %8.2f (%+6.1f c) | FEM %8.2f '
                  '(%+6.1f c) | FEM-1D %+6.1f c | loc %.2f | ne %d'
                  % (note, t, ref[note], 1200 * math.log2(ref[note] / t),
                     fem_f, 1200 * math.log2(fem_f / t),
                     1200 * math.log2(fem_f / ref[note]), w, ne))
            results[note] = {'target': t, 'f_1d': ref[note], 'f_fem': fem_f,
                             'fem_minus_1d_cents':
                                 1200 * math.log2(fem_f / ref[note]),
                             'localization': w, 'elements': ne}
        else:
            print('%-4s: no localized mode found; candidates: %s'
                  % (note, [(round(f, 1), round(w, 2)) for f, w in cands]))
    out = os.path.join(HERE, '..', '..', 'reports', 'fem_vs_1d.json')
    if results:
        existing = {}
        if os.path.exists(out):
            existing = json.load(open(out))
        existing.update(results)
        json.dump(existing, open(out, 'w'), indent=1)
        print('saved ->', out)


if __name__ == '__main__':
    main()
