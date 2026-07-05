# 3D Helmholtz gate for the M2d tapered whistle core (NGSolve/Netgen).
#
# Formulation: column-only with END-CORRECTION STUB terminations. Each open
# tonehole exit and the foot are extended by a cylindrical stub of length
# 0.6133*a (the unflanged end correction) terminated with p=0 — the same
# termination physics openwind's unflanged radiation applies in 1D, so
# FEM-vs-1D differences isolate the interior 3D effects under test:
# bore taper, graded crown depths, angled chimneys, hole-crown curvature.
# Closed holes end at the face plane as rigid walls (sealed).
#
# History: an exterior-air pocket harness was tried first and abandoned —
# four geometry-leak bugs later its cylindrical control never stabilized
# (see git history). The column-only formulation was clean in every probe.
#
# Usage:
#   .venv-fem/bin/python scripts/fem/m2d_helmholtz.py --note all
#   FEMGEO=ctl_fem_geometry.json ... --note all      (cylindrical control)
import argparse
import json
import math
import os

import numpy as np
from netgen.occ import Axes, Cone, Cylinder, Dir, OCCGeometry, Pnt, X
import ngsolve as ng

HERE = os.path.dirname(os.path.abspath(__file__))
GEONAME = os.environ.get('FEMGEO', 'm2d_fem_geometry.json')
G = json.load(open(os.path.join(HERE, '..', GEONAME)))

C_MM_S = 343_360.0
END_CORR = 0.6133           # unflanged end correction, x opening radius
HOLES = ['h1', 'h2', 'h3', 'h4', 'h5', 'h6']
BS = G['body_start']
LE = G['L_eff']


def r_bore(x):
    if x <= BS:
        return G['r_socket']
    f = (x - BS) / (LE - BS)
    return G['r_socket'] + f * (G['r_foot'] - G['r_socket'])


def face_plane():
    xs = [G['positions'][h] for h in ('h1', 'h6')]
    zs = [-(r_bore(G['positions'][h]) + G['crowns'][h]) for h in ('h1', 'h6')]
    slope = (zs[1] - zs[0]) / (xs[1] - xs[0])
    return zs[0] - slope * xs[0], slope


def path_to_plane(h):
    z0, sl = face_plane()
    hx = G['positions'][h]
    th = math.radians(G['theta_deg'][h])
    t = -(z0 + sl * hx) / (math.cos(th) + sl * math.sin(th))
    return t - r_bore(hx)


def build_domain(mask):
    """Bore + chimneys; open exits and the foot get end-correction stubs
    whose end discs (plus the window plane) are p=0."""
    nseg = 8
    xs = [0.0] + [BS + (LE - BS) * i / nseg for i in range(nseg + 1)]
    col = Cylinder(Pnt(0, 0, 0), X, r=G['r_socket'], h=BS)
    for i in range(nseg):
        x0, x1 = xs[i + 1], xs[i + 2]
        if abs(r_bore(x0) - r_bore(x1)) < 1e-9:
            col = col + Cylinder(Pnt(x0, 0, 0), X, r=r_bore(x0), h=x1 - x0)
        else:
            col = col + Cone(Axes(Pnt(x0, 0, 0), X, Dir(0, 0, 1)),
                             r_bore(x0), r_bore(x1), x1 - x0, 2 * math.pi)
    foot_ext = END_CORR * G['r_foot']
    col = col + Cylinder(Pnt(LE - 2, 0, 0), X, r=G['r_foot'], h=2.0 + foot_ext)
    ends = [(LE + foot_ext, 0.0, 'x')]
    for h in HOLES:
        hx = G['positions'][h]
        a = G['diameters'][h] / 2.0
        th = math.radians(G['theta_deg'][h])
        d = Dir(math.sin(th), 0, -math.cos(th))
        ln = path_to_plane(h) + r_bore(hx)
        if mask[h] == 'open':
            ln += END_CORR * a
            ends.append((hx + ln * math.sin(th), -ln * math.cos(th), 'p'))
        col = col + Cylinder(Pnt(hx, 0, 0), d, r=a, h=ln)
    for f in col.faces:
        c = f.center
        f.name = 'wall'
        if abs(c.x) < 1e-6:
            f.name = 'popen'
        for (ex, ez, kind) in ends:
            if kind == 'x' and abs(c.x - ex) < 0.05:
                f.name = 'popen'
            elif kind == 'p' and abs(c.x - ex) < 0.4 and abs(c.z - ez) < 0.4:
                f.name = 'popen'
    col.solids[0].maxh = 5.5
    return col


def eigs(mesh, fshift, nev=24):
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
        mesh = ng.Mesh(OCCGeometry(dom).GenerateMesh(maxh=5.5))
        mesh.Curve(2)
        fs = eigs(mesh, ref[n])
        t = G['targets'][n]
        near = [f for f in fs if abs(1200 * math.log2(f / ref[n])) < 100]
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
                  % (n, ref[n], [round(q, 1) for q in fs[:6]]), flush=True)
    out = os.path.join(HERE, '..', '..', 'reports',
                       'fem_gate_' + GEONAME.replace('.json', '') + '.json')
    json.dump(results, open(out, 'w'), indent=1)
    print('saved ->', out)


if __name__ == '__main__':
    main()
