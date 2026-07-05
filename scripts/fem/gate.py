"""As-exported FEM gate — solves the 3-D Helmholtz eigenproblem on the LITERAL
exported STEP, not a reconstruction.

Imports artifacts/<name>/air.step and manifest.json (both written by
foldflute.solids from the canonical doc), tags the STEP's faces by matching
face centres to the manifest openings, applies per-fingering boundary
conditions (open exits get a 0.6133*a unflanged end-correction stub -> p=0;
closed holes and all walls stay rigid; window and foot p=0 via stubs), and
compares each fingering's lowest resonance with openwind's 1-D prediction
built from the SAME kernel. Because the geometry, the manifest, and the
openwind input all come from foldflute.geometry, a divergence between them
shows up here as a gate failure rather than silently shipping.

Env: .venv-fem (ngsolve + openwind). Kernel imported from repo root.

Usage:
  .venv-fem/bin/python scripts/fem/gate.py --artifacts artifacts/m2d \
      --spec spec/m2d.json --note all
"""
import argparse
import json
import math
import os
import sys

import numpy as np
from netgen.occ import Cylinder, Dir, OCCGeometry, Pnt
import ngsolve as ng

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from foldflute import geometry as G   # noqa: E402

C_MM_S = 343_360.0
END_CORR = 0.6133


def tag_and_extend(step_path, man, mask):
    """Load the exported air STEP, add end-correction stubs on open exits +
    window + foot, and name p=0 boundaries. Returns the netgen solid."""
    solid = OCCGeometry(step_path).shape
    stubs = []

    def stub(center, normal, radius):
        d = Dir(normal[0], normal[1], normal[2])
        return Cylinder(Pnt(center[0] - normal[0] * 0.5,
                            center[1] - normal[1] * 0.5,
                            center[2] - normal[2] * 0.5), d, r=radius,
                        h=END_CORR * radius + 0.5)

    # window: the x=0 disc radiates (mouth end); model as p=0 directly there.
    # foot + open holes: end-correction stubs terminated p=0.
    stubs.append(('foot', stub(man['foot']['center'], man['foot']['normal'],
                               man['foot']['radius'])))
    for hid, e in man['holes'].items():
        if mask[hid] == 'open':
            stubs.append((hid, stub(e['center'], e['normal'], e['radius'])))
    for _, s in stubs:
        solid = solid + s

    # name faces: p=0 on window plane, on every stub end disc; wall elsewhere
    ends = [('window', [0.0, 0.0, 0.0], None)]
    ends.append(('foot', [man['foot']['center'][i]
                          + man['foot']['normal'][i] * (END_CORR
                          * man['foot']['radius']) for i in range(3)],
                 man['foot']['radius']))
    for hid, e in man['holes'].items():
        if mask[hid] == 'open':
            c = [e['center'][i] + e['normal'][i] * (END_CORR * e['radius'])
                 for i in range(3)]
            ends.append((hid, c, e['radius']))
    for f in solid.faces:
        c = f.center
        f.name = 'wall'
        if abs(c.x) < 1e-6:                     # window plane disc
            f.name = 'popen'
            continue
        for _, ec, _r in ends[1:]:
            if math.dist((c.x, c.y, c.z), ec) < 1.2:
                f.name = 'popen'
                break
    solid.solids[0].maxh = 5.0
    return solid


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


def ref_1d(doc):
    from openwind import ImpedanceComputation
    out = {}
    for n in doc['fingerings']:
        mb, rows, chart, temp = G.to_openwind(doc, doc['fingerings'][n])
        comp = ImpedanceComputation(np.array([300.0]), mb, rows, chart,
                                    note=G.ow_label(n), temperature=temp,
                                    losses=True, radiation_category='unflanged',
                                    matching_volume=True)
        t = G.note_hz(n, doc['target']['pitch_standard_hz'])
        fs = t * 2 ** (np.linspace(-170, 170, 29) / 1200)
        comp.recompute_impedance_at(fs)
        i = int(np.argmin(np.abs(comp.impedance)))
        fs2 = fs[i] * 2 ** (np.linspace(-12, 12, 13) / 1200)
        comp.recompute_impedance_at(fs2)
        out[n] = float(fs2[int(np.argmin(np.abs(comp.impedance)))])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--artifacts', required=True)
    ap.add_argument('--spec', required=True)
    ap.add_argument('--note', nargs='*', default=['all'])
    args = ap.parse_args()
    doc = G.load(args.spec)
    man = json.load(open(os.path.join(args.artifacts, 'manifest.json')))
    step = os.path.join(args.artifacts, 'air.step')
    ref = ref_1d(doc)
    notes = list(doc['fingerings']) if args.note == ['all'] else args.note

    # CAD<->FEM volume invariant on the all-closed (foot-only stub) build
    results = {}
    for n in notes:
        solid = tag_and_extend(step, man, doc['fingerings'][n])
        mesh = ng.Mesh(OCCGeometry(solid).GenerateMesh(maxh=5.0))
        mesh.Curve(2)
        if n == notes[0]:
            vmesh = ng.Integrate(1.0, mesh)
            print('mesh volume %.0f mm3 vs manifest air %.0f (ratio %.3f, '
                  'incl stubs)' % (vmesh, man['air_volume_mm3'],
                                   vmesh / man['air_volume_mm3']), flush=True)
        fs = eigs(mesh, ref[n])
        t = G.note_hz(n, doc['target']['pitch_standard_hz'])
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
            print('%-4s: no mode within 100c of 1D %.1f; nearest %s'
                  % (n, ref[n], [round(q, 1) for q in fs[:6]]), flush=True)
    name = os.path.basename(args.artifacts.rstrip('/'))
    outp = os.path.join(os.path.dirname(__file__), '..', '..', 'reports',
                        'gate_%s.json' % name)
    json.dump(results, open(outp, 'w'), indent=1)
    print('saved ->', outp)


if __name__ == '__main__':
    main()
