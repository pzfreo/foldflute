"""Regenerate all M2d artifacts from the canonical doc and run invariants.

One command turns spec/m2d.json into artifacts/m2d/{air.step, body.step,
manifest.json} and asserts the generator did not silently produce wrong
geometry. Fails loud (non-zero exit) on any invariant breach. The FEM gate
(scripts/fem/gate.py) is the separate, heavier acoustic check run afterward.

Env: .venv-cad.  Usage: .venv-cad/bin/python scripts/regen.py [spec/m2d.json]
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from build123d import import_step                       # noqa: E402
from foldflute import geometry as G, solids as S        # noqa: E402


def main():
    spec = sys.argv[1] if len(sys.argv) > 1 else 'spec/m2d.json'
    out = os.path.join('artifacts', G.load(spec)['name'])
    doc = G.load(spec)                                   # re-validates
    fails = []

    air_step, air_vol = S.export_air(doc, out)
    body_step, body_vol = S.export_body(doc, out)
    print('generated %s (air %.0f mm3) and %s (body %.0f mm3)'
          % (air_step, air_vol, body_step, body_vol))

    # 1. STEP round-trip: the shipped artifact must re-import identically
    ra = import_step(air_step).volume
    if abs(ra / air_vol - 1.0) > 1e-4:
        fails.append('air.step round-trip volume drift %.5f' % (ra / air_vol))

    # 2. analytic tripwire: catch gross generator failure (crude ±12%)
    an = G.analytic_air_volume(doc)
    if abs(air_vol / an - 1.0) > 0.08:
        fails.append('air volume %.0f vs analytic %.0f off >8%%'
                     % (air_vol, an))
    print('  round-trip match %.5f | analytic ratio %.3f'
          % (ra / air_vol, air_vol / an))

    # 3. within-hand gap constraint (spec v0.4 hard limit)
    limit = doc['constraint']['max_within_hand_gap_mm']
    gaps = G.within_hand_gaps(doc)
    for pair, g in gaps.items():
        if g > limit + 1e-6:
            fails.append('within-hand gap %s = %.1f > %.1f mm'
                         % (pair, g, limit))
    print('  within-hand gaps ' + '  '.join('%s %.1f' % (k, v)
                                             for k, v in gaps.items())
          + '  (limit %.0f)' % limit)

    # 4. body must be a single closed solid (printability precondition)
    b = import_step(body_step)
    nsolids = len(b.solids()) if hasattr(b, 'solids') else 1
    if nsolids != 1:
        fails.append('body is %d solids, expected 1' % nsolids)

    if fails:
        print('\nINVARIANTS FAILED:')
        for f in fails:
            print('  x ' + f)
        sys.exit(1)
    print('\nall invariants PASS')


if __name__ == '__main__':
    main()
