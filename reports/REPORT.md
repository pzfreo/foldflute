# foldflute first end-to-end run — 2026-07-03

Toolchain run against SPEC.md draft v0.3, driven through `build123d-mcp`
(0.3.62, augura 0.1.5, openwind via `uv tool run`). Deliverables: solved M1
straight and M2 folded designs, gated print files in `artifacts/`, solver
sources in `foldflute/prelude.py`, parametric build scripts in `scripts/`.

## 1. Acoustic model

Openwind 1D FEM, viscothermal losses on, unflanged radiation at holes and
foot, 20 °C, impedance **minima** at the window datum (x=0), headjoint
effective-length correction 0.0 (to be calibrated in M1). Resonances located
by a local geometric frequency grid ±170 cents around each target with
two-stage parabolic refinement — openwind's phase-based
`antiresonance_frequencies()` was unreliable on coarse grids.

Hole positions and bore length are solved by Gauss–Seidel (bore end → D4,
then h6→h1 for E4…C#5) with damped secant steps; register 1 converges to
±0.1 cents in 2–3 sweeps. Register 2 and D6 are *predictions*, controlled by
the diameter choice, which is an input (see §2).

## 2. Diameter scenario study (straight bore)

| scenario | diameters h1..h6 | worst R2 note | h1–h6 span | fold-layout gain |
|---|---|---|---|---|
| S0 seed | 8/8/8/8/8.5/9 | C#6 −97 c | 209 mm | — |
| **S1 big (selected)** | 11/10.5/10/10/10.5/11 | C#6 −26 c | 205 mm | none |
| S2 compression | 11.5/11/10.5/8/7/6.5 | F#5 −37 c | 180 mm | ~25 mm |
| S3 compromise | 11.5/11/10.5/9.5/9/8.5 | C#6 −24 c | 194 mm | ~10 mm |

Selected S1: large holes buy ~70 cents on the top of register 2 (the
open-lattice cutoff moves up); diameter "compression" of the hole span costs
real register-2 accuracy and recovers only ~25 mm of the ~110 mm needed to
reach the low A layout. Passive-model caveat: players overblow the second
register slightly sharp, so the −10…−26 c predictions up high should play
closer than they read.

## 3. Solved designs

Both use S1 diameters, bore Ø22 cylindrical, chimneys 3.30 mm (measured from
the as-built CAD wall and written back per spec; the 3.0 mm estimate shifted
holes by ~1 mm, above the rebuild threshold, so geometry was rebuilt).

**M1 straight**: bore end 571.60 mm from window. Holes (mm from window):
283.6 / 319.3 / 358.9 / 411.1 / 431.9 / 488.3.

**M2 folded**: bore end 574.57 (+2.97 vs M1 — the bend correction re-solve),
bend start 464.90, h6 at 41.8° into the bend. Holes: 283.6 / 319.3 / 358.9 /
411.1 / 431.9 / 489.0.

Model cents (register 1 all ±0.4 c by construction):

| note | D4 | E4 | F#4 | G4 | A4 | B4 | C#5 | D5 | E5 | F#5 | G5 | A5 | B5 | C#6 | D6 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| M1 | 0.0 | −0.1 | −0.2 | −0.1 | +0.4 | −0.4 | −0.1 | +7.7 | −2.1 | −9.5 | −6.7 | −17.7 | −23.0 | −25.8 | +13.1 |
| M2−M1 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |

M2 matches M1 within 0.03 cents per fingering in-model (the M2 exit
criterion of ±10 c is then a test of the bend-correction model itself once
printed). G3's ±15 c pre-calibration target is met for register 1 + D5/E5/
F#5/G5/D6; A5/B5/C#6 sit at −18…−26 c passive.

## 4. Bend correction (SPEC 6.2)

Toroidal bend, long-wavelength limit, potential (1/r) flow across the
section: inertance ratio bend/straight = (1+√(1−B²))/2 ≈ 1−B²/4, B = bore
radius / bend radius. Same volume ⇒ compliance unchanged ⇒ the bend is
acoustically shorter than its arc by that fraction (pitch rises). Applied as
a piecewise arc-length remapping with κ=1 (velocity-antinode worst case;
κ is the calibration knob). For R=33/bore 11: 2.83% of the 103.7 mm bend arc
= 2.9 mm, ≈ +9 c on D4 — matching the observed +2.97 mm bore-end re-solve.
Refs: Nederveen JASA 104(3):1616 (1998); Félix, Dalmont, Nederveen JASA
131:4164 (2012). (Exact published coefficients are paywalled; the derivation
above is in `foldflute/prelude.py` and is calibration-refinable.)

## 5. Fold-path result — the honest news

With a Ø22 cylindrical bore, a 180° bend needs πR ≥ 104 mm of arc at the
minimum legal radius (1.5×bore Ø = 33 mm):

- between adjacent toneholes: < 57 mm available → impossible;
- above h1 (order-preserving fold-back = 2 bends ≥ 207 mm): ~80 mm of
  straight available below the 56 mm spigot → impossible;
- wholly below h6: ~90 mm of tail → impossible.

**The only legal fold is a foot U-bend with h6 partway into the bend**
(41.8°, 10.5 mm lateral offset — comparable to offset pinky holes on
flutes). Solved fold: R=33, passes 66 mm apart in-plane, visible body 338.5
mm vs low A target 286.4 (+52.1, relaxed per `relax_visual_length`); foot
opening returns up the second pass 6 mm from the bend apex (the low point —
condensate drains to the opening; no drain hole needed).

The low A hole layout itself is **infeasible for any fold topology**: hole
spacing on a single pass equals acoustic spacing (span ~205 mm vs low A's
139.5 mm), and diameter compression tops out ~25 mm (§2). Exterior errors vs
the MK Pro A targets: h1 +17, h2 +27, h3 +40, h4 +58, h5 +59, h6 +81 mm
(tolerance 8 mm). Per SPEC §7 the layout is relaxed rather than hand-editing
acoustics. Datum assumption: `hole_centers_from_socket` measured from where
the body emerges from the sleeve (arc 176); if the datum is the tube end
(arc 120), all errors grow by 56 mm — worth confirming against the MK
drawing.

A 4-pass "S above the holes + foot U" topology could hit the 286 mm length
but needs ≥207 mm folded above h1 where only ~80 mm exists, and would push
h1 above the socket. Options if the low A length matters more than v0.3's
constraints allow: taper/neck the foot bore (relaxes local bend-radius
floor), allow smaller top-fold radii, or accept a ~95 mm-wide 3-tube bundle
with bassoon-style hole ergonomics. Parked for discussion.

## 6. Print architecture and gates

| part | geometry | print | Augura |
|---|---|---|---|
| M2 section_a | x 0–172 + 20 tenon, carries h1, spigot Ø24.8 | vertical, monolithic | 4 warnings¹ |
| M2 left_shell (front) | +z half of section B; carries h2–h6 wholly (seam never crosses a hole) | seam face down | 2 warnings + 3 info |
| M2 right_shell (back) | −z half | seam face down | same |
| M2 alignment_pin ×3 | Ø2.3×4 loose dowels into Ø2.5 recesses in 3 seam bosses | standing | brim advisory |
| M1 section_a | x 0–219 + tenon, h1–h2 | vertical | 5 warnings¹ |
| M1 section_b | x 219–451.6, h3–h6, socket ring | vertical | 3 warnings¹ |

¹ All findings are warnings/info, no failures: (a) tenon/socket stop-ring
ceilings (1.3–2.1 mm ring bridges — left as-is, they are the insertion
stops); (b) tonehole ceilings on vertical prints (4.6–4.8 mm bridges —
routine, ream after printing if sagged); (c) brim recommended on tall
sections (aspect 15–19 — yes, use a brim); (d) 0.2–0.4 mm "thin wall" =
tangential feather edges where hole cylinders exit curved walls, inherent to
drilled toneholes. Augura reported no overhang findings on the Ø22 bore
crowns of the flat-printed shells; treat crown sag as a print-quality watch
item regardless.

Section A is deliberately monolithic (deviation from a whole-body shell
split): the two sealing surfaces — headjoint spigot and joint tenon — stay
seamless and round. Joints: tenon Ø26.0 into socket Ø26.2 (0.2 diametral,
same as the headjoint sleeve fit), 20 mm engagement, shoulder stop, waxed
thread fallback per spec.

All fit checks passed: tenon/socket touching with zero interpenetration,
shells touching on the full seam, air column touching walls with zero
intersection volume. All parts pass the CAD validity gate (watertight
manifold single solids). Footprints ≤ 227×101 mm (bed 256², limit 250²).

## 7. Session evidence

- `reports/solved_designs.json` — full solved specs + cents tables + fold report
- `reports/solver_run_summary.json` — scenario study + provenance
- `reports/m2_*.png` — renders (iso, assembly, clipped)
- `scripts/session_build_log.py` — raw executed-session log (provenance)
- `scripts/build_m*.py` — clean parametric rebuilds
- `artifacts/m1/`, `artifacts/m2/` — STLs (bed-oriented) + canonical STEPs

## 8a. Addendum — convoluted internal bore study (same day)

Clarified intent: the goal was never "shorten with a bend" but a **convoluted
internal bore** so the exterior compresses to low A finger positions. This was
studied quantitatively (`reports/convolution_feasibility.json`, solver code in
`foldflute/prelude.py` v2 section). Result: **infeasible at a Ø22 cylindrical
bore with plain keyless toneholes**, from three directions:

1. **Geometry, shallow chimneys.** Hole span is pitch-set (~205 mm for S1
   diameters; 180 mm at the acoustically-damaged S2 extreme) vs low A's
   139.5 mm. Toroidal detours need ~4R ≈ 132 mm of axial footprint per
   inter-hole gap of 20–34 mm — impossible at any bore radius. Mitred folds
   (Coltman 2006: compensated beveled miter = straight tube −0.32×ID per 90°,
   mode-independent — adopted as the correction model for any internal fold)
   fit *outside* the hole span but cannot compress spacing *within* it.
2. **Stepped deep chimneys.** Holes tapping runs at different depths shift by
   measured −70…−118 cents (−30…−55 mm) — adjacent 25 mm-spaced holes invert
   their physical order. Dead end.
3. **Graded deep chimneys ("ramp bore").** A smooth chimney ramp (10→43 mm,
   no elbows in the hole span) genuinely places all six holes within
   **9 mm of the low A targets** with register 1 in tune — but register 2
   separates catastrophically (F#5 −147 c, G5 −158 c, A5 +128 c): each deep
   tonehole becomes a resonant side branch and the instrument cannot
   overblow. A mild ramp (4→18 mm) is dominated: span only 175 mm, worst
   layout error still +50 mm, register 2 already at −75 c.

The trade frontier (span vs worst register-2 error): 205 mm/−26 c (flat),
175 mm/−75 c (mild ramp), 155 mm/−158 c (full ramp). Register 2 degrades
faster than the layout improves at every point.

**Paths that could still reach low A hand feel** (owner's call, all spec
changes):

- **(a) Ship v1 as-is** — foot-fold M2, proper two-register whistle, length
  339 vs 286 target, hole positions low-D-like.
- **(b) Low A *length*, shifted hand.** Fold ~58 mm of bore into a deep
  gallery above h1 (one mitred sidestep pair fits the 46 mm footprint) plus
  the tail gallery below h6, keeping all chimneys shallow: full 286 mm low A
  body length, S1 register quality, hole block sits ~40 mm higher than the
  low A hand with low-D spacing. Feasible with today's solvers; not yet
  CAD-built.
- **(c) Keys/touch-plates** at low A finger positions actuating pads over
  acoustically-placed holes — the way every large woodwind solves exactly
  this. Printable rocker keys are plausible FDM parts; significant spec/scope
  change (v0.3 non-goal).
- **(d) Narrower/tapered bore study** — second-order for span (hole spacing
  is pitch-set), changes the instrument's voice; not recommended as the
  primary lever.

## 8. Open items / assumptions to confirm

1. **Headjoint effective length 0.0** — M1 exists to calibrate this; expect
   all notes to shift together on the first print (fit with
   `global_length_correction` + `headjoint_effective_length_correction`).
2. Low A datum ambiguity (§5).
3. Insertion depth assumed full (56 mm): body acoustic start = tube end.
4. Register-2 passive flatness (A5/B5/C#6): revisit after M1 tuner data;
   options are bigger h1/h2, chimney thinning at the top holes, or accepting
   player compensation.
5. `H6_MAX_BEND_ANGLE_DEG = 50°` is an ergonomic judgment — verify R3 can
   seal a hole offset 10.5 mm laterally on the printed prototype.
6. Bend-correction κ=1.0 — M2-vs-M1 tuner comparison calibrates it.
7. `apply_calibration()` (fitting corrections from tuner observations) is
   not yet implemented — deliberately deferred until M1 data exists.
