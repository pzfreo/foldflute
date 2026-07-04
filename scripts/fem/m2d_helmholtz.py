# 3D Helmholtz re-gate of the M2d tapered-serpentine acoustic core.
#
# Gates what is NEW versus the already-FEM-validated cylindrical build:
# the 22->14mm bore taper, the curved bore line (graded crown depth), and
# the deep angled chimneys. The bore axis is modelled unrolled (straight):
# the fold/gallery elbows were validated separately (455 vs 452mm effective).
# Bore = piecewise-conical frustums; face = single tilted plane through the
# chimney exits; each chimney is a tilted cylinder of the solved angle.
# Exterior air block gives real radiation loading; p=0 at window plane and
# far boundary; fingertip caps seal closed holes.
#
# Usage: .venv-fem/bin/python scripts/fem/m2d_helmholtz.py --note all
import argparse
import json
import math
import os

import numpy as np
from netgen.occ import (Axes, Box, Cone, Cylinder, Dir, Glue,
                        OCCGeometry, Pnt, X, Z)
import ngsolve as ng

HERE = os.path.dirname(os.path.abspath(__file__))
G = json.load(open(os.path.join(HERE, '..', 'm2d_fem_geometry.json')))

C_MM_S = 343_360.0
HOLES = ['h1', 'h2', 'h3', 'h4', 'h5', 'h6']
BS = G['body_start']          # 120: taper begins here; window at x=0
LE = G['L_eff']


def r_bore(x):
    if x <= BS:
        return G['r_socket']
    f = (x - BS) / (LE - BS)
    return G['r_socket'] + f * (G['r_foot'] - G['r_socket'])


def face_plane():
    """Tilted top face through the chimney exit points (z negative side).
    Returns (z0, slope): z_face(x) = z0 + slope*x."""
    xs = [G['positions'][h] for h in ('h1', 'h6')]
    zs = [-(r_bore(G['positions'][h]) + G['crowns'][h]) for h in ('h1', 'h6')]
    slope = (zs[1] - zs[0]) / (xs[1] - xs[0])
    return zs[0] - slope * xs[0], slope


def path_to_plane(h):
    """Chimney path length from bore wall to the tilted face plane."""
    z0, sl = face_plane()
    hx = G['positions'][h]
    th = math.radians(G['theta_deg'][h])
    t = -(z0 + sl * hx) / (math.cos(th) + sl * math.sin(th))
    return t - r_bore(hx)


def air_column():
    """Piecewise-conical bore (axis on z=0) + angled chimney tubes."""
    nseg = 8
    xs = [0.0] + [BS + (LE - BS) * i / nseg for i in range(nseg + 1)]
    col = Cylinder(Pnt(0, 0, 0), X, r=G['r_socket'], h=BS)
    for i in range(nseg):
        x0, x1 = xs[i + 1], xs[i + 2]
        col = col + Cone(Axes(Pnt(x0, 0, 0), X, Dir(0, 0, 1)), r_bore(x0), r_bore(x1),
                         x1 - x0, 2 * math.pi)
    z0, sl = face_plane()
    for h in HOLES:
        hx = G['positions'][h]
        th = math.radians(G['theta_deg'][h])
        d = Dir(math.sin(th), 0, -math.cos(th))   # toward the face (z<0)
        ln = path_to_plane(h) + r_bore(hx) + 5.0
        col = col + Cylinder(Pnt(hx, 0, 0), d, r=G['diameters'][h] / 2.0, h=ln)
    # trim 0.6mm BELOW the face plane: chimney stubs protrude into the
    # exterior so the col/ext boolean imprints real shared faces
    trim = lower_halfspace(z0 - 0.6, sl) * Box(
        Pnt(G['positions']['h1'] - 25, -250, -500), Pnt(1000, 250, 100))
    col = col - trim
    # foot stub: 1mm past the body end face
    col = col + Cylinder(Pnt(LE - 5, 0, 0), X, r=G['r_foot'], h=6.0)
    return col


def lower_halfspace(z0, sl):
    """Solid occupying the region below the tilted face plane z=z0+sl*x."""
    from netgen.occ import Axis as OAxis
    b = Box(Pnt(-400, -250, -500), Pnt(1000, 250, 0))
    b = b.Rotate(OAxis(Pnt(0, 0, 0), Dir(0, 1, 0)),
                 -math.degrees(math.atan(sl)))
    return b.Move((0, 0, z0))


def build_domain(mask):
    z0, sl = face_plane()
    col = air_column()
    # body: slab around the bore, top face = tilted plane
    body = Box(Pnt(BS - 60, -15, -40), Pnt(LE, 15, 25))
    body = body - lower_halfspace(z0, sl)
    ext = Box(Pnt(BS - 40, -80, -95), Pnt(LE + 60, 80, 60))
    ext = ext - body
    for h in HOLES:
        if mask[h] == 'closed':
            hx = G['positions'][h]
            th = math.radians(G['theta_deg'][h])
            d = Dir(math.sin(th), 0, -math.cos(th))
            ln0 = path_to_plane(h) + r_bore(hx)
            ext = ext - Cylinder(Pnt(hx, 0, 0), d,
                                 r=G['diameters'][h] / 2.0 + 4.0,
                                 h=ln0 + 5.0)
    ext = ext - col
    col.solids[0].name = 'col'
    col.solids[0].maxh = 5.5
    for s in ext.solids:
        s.name = 'ext'
        s.maxh = 16.0
    dom = Glue([col, ext])
    for f in dom.faces:
        c = f.center
        if abs(c.x) < 1e-6:
            f.name = 'popen'
        elif (abs(c.x - (BS - 40)) < 1e-6 or abs(c.x - (LE + 60)) < 1e-6
              or abs(abs(c.y) - 80) < 1e-6 or abs(c.z + 95) < 1e-6
              or abs(c.z - 60) < 1e-6):
            f.name = 'popen'
        else:
            f.name = 'wall'
    return dom


def eigs(mesh, fshift, nev=140):
    fes = ng.H1(mesh, order=2, dirichlet='popen')
    u, v = fes.TnT()
    a = ng.BilinearForm(ng.grad(u) * ng.grad(v) * ng.dx).Assemble()
    m = ng.BilinearForm(u * v * ng.dx).Assemble()
    k0 = 2 * math.pi * fshift / C_MM_S
    vecs = [a.mat.CreateColVector() for _ in range(nev)]
    lams = ng.ArnoldiSolver(a.mat, m.mat, fes.FreeDofs(), vecs, shift=k0 * k0)
    return sorted(math.sqrt(l.real) * C_MM_S / (2 * math.pi)
                  for l in lams if l.real > 0)


def ref_1d():
    from openwind import ImpedanceComputation
    out = {}
    for n, fing in G['fingerings'].items():
        main = [[0.0, G['r_socket'] / 1000], [BS / 1000, G['r_socket'] / 1000],
                [LE / 1000, G['r_foot'] / 1000]]
        rows = [['label', 'position', 'radius', 'chimney']]
        for h in HOLES:
            t = path_to_plane(h) + (G['dt_open']
                                    if fing[h] == 'open' else 0.0)
            rows.append([h, G['positions'][h] / 1000,
                         G['diameters'][h] / 2000, t / 1000])
        labels = [x.replace('#', 's') for x in G['fingerings']]
        chart = [['label'] + labels]
        for h in HOLES:
            chart.append([h] + ['x' if G['fingerings'][m][h] == 'closed'
                                else 'o' for m in G['fingerings']])
        comp = ImpedanceComputation(np.array([300.0]), main, rows, chart,
                                    note=n.replace('#', 's'),
                                    temperature=20.0, losses=True,
                                    radiation_category='unflanged',
                                    matching_volume=True)
        t = G['targets'][n]
        fs = t * 2 ** (np.linspace(-170, 170, 29) / 1200)
        comp.recompute_impedance_at(fs)
        i = int(np.argmin(np.abs(comp.impedance)))
        fs2 = fs[i] * 2 ** (np.linspace(-12, 12, 13) / 1200)
        comp.recompute_impedance_at(fs2)
        out[n] = float(fs2[int(np.argmin(np.abs(comp.impedance)))])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--note', nargs='*', default=['all'])
    args = ap.parse_args()
    ref = ref_1d()
    notes = list(G['fingerings']) if args.note == ['all'] else args.note
    results = {}
    for n in notes:
        dom = build_domain(G['fingerings'][n])
        mesh = ng.Mesh(OCCGeometry(dom).GenerateMesh(maxh=16.0))
        mesh.Curve(2)
        fs = eigs(mesh, ref[n])
        t = G['targets'][n]
        near = [f for f in fs
                if abs(1200 * math.log2(f / ref[n])) < 100]
        if near:
            f = min(near, key=lambda q: abs(q - ref[n]))
            print('%-4s target %8.2f | 1D %8.2f (%+6.1f c) | FEM %8.2f '
                  '(%+6.1f c) | FEM-1D %+6.1f c | ne %d'
                  % (n, t, ref[n], 1200 * math.log2(ref[n] / t), f,
                     1200 * math.log2(f / t), 1200 * math.log2(f / ref[n]),
                     mesh.ne), flush=True)
            results[n] = {'target': t, 'f_1d': ref[n], 'f_fem': f,
                          'fem_minus_1d_cents': 1200 * math.log2(f / ref[n])}
        else:
            print('%-4s: no mode within 100c of 1D %.1f; nearest: %s'
                  % (n, ref[n], [round(q, 1) for q in fs[:8]]), flush=True)
    out = os.path.join(HERE, '..', '..', 'reports', 'm2d_fem_vs_1d.json')
    json.dump(results, open(out, 'w'), indent=1)
    print('saved ->', out)


if __name__ == '__main__':
    main()
