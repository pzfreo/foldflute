# foldflute geometry pipeline

One document, one kernel, one CAD path, one gate. Built to remove the failure
mode where geometry was described several times (solver state, ad-hoc JSON,
hand-written FEM reconstruction, hand-built preview) and the acoustic gate ran
on a *reconstruction* rather than the shipped artifact.

## Data flow

```
spec/m2d.json                         canonical document (edit this)
      │
      ▼
foldflute/geometry.py   pure kernel (math+json; imports in every venv)
      │  r_bore · face_plane · hole_exit · manifest · to_openwind · invariants
      ├────────────────────────────┐
      ▼                            ▼
foldflute/solids.py          scripts/fem/gate.py
 (.venv-cad, build123d)        (.venv-fem, ngsolve)
   air.step + body.step          imports air.step + manifest.json,
   + manifest.json               tags faces, 3-D Helmholtz per fingering,
      │                          compares to openwind 1-D from the SAME kernel
      ▼
scripts/regen.py  →  artifacts/m2d/  + invariant checks (fails loud)
```

The kernel is the *only* place geometric numbers are computed. `solids.py` and
`gate.py` both derive from it, so CAD, 1-D acoustics, and 3-D FEM cannot
silently describe different geometry. The FEM meshes the **literal exported
STEP** (not a mirror); `manifest.json`, also kernel-derived, tells it which
face is the window / each hole exit / the foot.

## Commands

```sh
# regenerate all artifacts from the doc + run invariants (build123d venv)
.venv-cad/bin/python scripts/regen.py                # -> artifacts/m2d/, PASS/FAIL

# as-exported acoustic gate on the shipped STEP (ngsolve venv)
.venv-fem/bin/python scripts/fem/gate.py \
    --artifacts artifacts/m2d --spec spec/m2d.json --note all
```

## Invariants (regen, all must PASS before print)

1. **round-trip** — air.step re-imports at volume match 1.00000 (export sane)
2. **analytic tripwire** — CAD air volume within 8% of the kernel's analytic
   integral (catches gross generator failure; currently 0.994)
3. **within-hand gaps** — every constrained hole pair ≤ 29 mm c-c (spec v0.4)
4. **single solid** — body.step is one closed solid (printability precondition)

The gate adds the cross-environment check: **mesh volume vs manifest air
volume** (currently 1.004 — same solid on both sides), then the acoustic
FEM-vs-1D table per fingering.

## Environments

- `.venv-cad` — build123d + openwind + scipy + numpy (solids, 1-D acoustics,
  the solver re-loop)
- `.venv-fem` — ngsolve + openwind + numpy (the gate)
- both import `foldflute/geometry.py` unchanged; the build123d-mcp session is
  now only an interactive inspector (import the STEP to render/measure), not a
  build authority.

## Status / next

The doc currently holds the pre-correction M2d solve. The gate confirms the
pipeline reproduces the earlier reconstruction's differentials (up-bore-tilted
h3/h6 read sharp; down-bore tilts clean; ~+24 c all-closed baseline offset).
Next design step: fit per-hole corrections from `reports/gate_m2d.json`,
re-solve into `spec/m2d.json`, `regen`, gate again — one loop, one truth. Then
the print-part CAD (sections/joint/chamfers) extends `solids.py`.
Known rough edge: the gate's mode picker occasionally mis-selects on all-open
fingerings (C#5); needs a mode-shape (energy-localization) filter.
