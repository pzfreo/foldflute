# foldflute - Folded Low-D Whistle Designer

**Project specification - draft v0.3 - 2026-07-03**  
**Owner:** Paul Fremantle

## 1. Overview

foldflute is a single-user desktop toolchain for designing a **folded-bore low D
tin whistle body** that:

- sounds as a low D whistle, with D4 as the fundamental;
- uses **standard six-hole tin whistle fingerings**;
- reuses an existing, known-good whistle headjoint for the windway, window, blade,
  and voicing-critical geometry;
- has the **physical length and approximate finger-hole layout of a low A
  whistle body** supplied by the user.

The bore is folded so that acoustic hole positions along the air column can land
near low A-style exterior finger positions. The acoustic target remains low D;
"low A body size" means physical length and hand layout, not a low A pitch
target.

An LLM orchestrates the design loop through `build123d-mcp` as the primary CAD
workbench: it builds geometry incrementally, renders views, measures dimensions,
checks clearances, compares shape revisions, exports artifacts, and runs
printability analysis. Deterministic solvers own the acoustically critical
numbers: bore length, tonehole positions, tonehole diameters, bend corrections,
chimney lengths, and printability checks.

**Core architectural insight:** the passive pitch model depends on the 1D air
column: bore radius vs. acoustic arc length, side-hole geometry, radiation
conditions, and boundary corrections. The exterior 3D body shape is mostly
acoustically free, provided the CAD-derived chimney depths and bend corrections
are fed back into the acoustic model.

## 1a. Requirement amendment — v0.4 (owner, 2026-07-05)

The low A physical layout goals are **withdrawn**. `physical_layout`
hole-position targets and the low A body-length target no longer constrain
the design. The binding ergonomic requirement is:

- **Adjacent holes within each hand are <= 29.0 mm centre to centre**
  (h1-h2, h2-h3 for the left hand; h4-h5, h5-h6 for the right hand).
- The fipple-to-left-hand distance and the left-hand-to-right-hand distance
  are unconstrained.

Consequences: G4 and the fold-path/layout machinery become optional
(compactness may still motivate folding, but nothing requires it); the
between-hands gap floats to whatever acoustics prefers; body length floats.

## 2. Goals

- G1. Produce a low D whistle body for an existing fixed headjoint.
- G2. Preserve standard six-hole tin whistle fingerings for the D major scale.
- G3. Predict resonance tuning before printing, accurate enough that one
  calibration print should suffice: relative tuning within +/-15 cents before
  calibration and +/-10 cents after one calibration iteration.
- G4. Fold the bore so acoustic hole positions map as closely as possible to a
  user-supplied low A physical layout: overall body length and approximate
  exterior hole/finger positions.
- G5. Wrap the air column in a printable body, with verified wall thickness,
  clearances, drainage, split-shell features, Augura printability analysis, and
  STL export.
- G6. Keep the full loop drivable by an LLM through `build123d-mcp`, using MCP
  tools for normal CAD iteration as well as final verification, while enforcing
  invariants in deterministic helper code.

## 3. Non-goals

- **Headjoint or voicing design.** Windway height, window length, blade geometry,
  chamfers, ramp shape, and mouthpiece ergonomics are copied from the existing
  headjoint and are not optimized.
- **Changing whistle fingering semantics.** v1 targets the standard D whistle
  six-hole diatonic chart. Chromatic half-holing and alternate cross-fingerings
  may be recorded later, but they are not optimized in v1.
- **Transverse flutes or reed instruments.** This project is for fipple whistles.
- **Sound synthesis or playability simulation.** We predict passive resonances,
  not tone, chiff, breath curve, or response.
- **Hosted, multi-user, or cloud service concerns.** This is a local desktop
  workflow.

## 4. System Architecture

```
                 +------------------------------+
                 | InstrumentSpec (JSONC)       |  single source of truth
                 +------+-----------+-----------+
                        |           |
          unrolled 1D   |           |  folded centerline + body params
                        v           v
      +----------------------+   +--------------------------+
      | Acoustic engine       |   | Geometry engine           |
      | (openwind helpers in  |   | (build123d via MCP)       |
      |  the CAD session)     |   |                          |
      | - impedance minima    |   | - folded bore sweep       |
      | - cents per fingering |   | - headjoint interface     |
      | - hole optimizer      |   | - tonehole chimneys       |
      | - bend corrections    |   | - shell split/export      |
      +----------+-----------+   | - clearance checks        |
                 |               | - Augura printability     |
                 |               +------------+-------------+
                 | measured chimney depths                 |
                 v                                         v
      +----------------------+   +--------------------------+
      | Fold-path solver      |   | Calibration protocol      |
      | (scipy constraints    |   | (print -> tuner -> fit    |
      |  + low A layout fit)  |   |  correction -> reprint)   |
      +----------+-----------+   +--------------------------+
                 ^
                 |
        +--------+---------+
        | LLM + human      | proposes topology/body shape,
        +------------------+ reads cents tables, drives iteration
```

### Division of Labour

| Concern | Owner |
|---|---|
| Headjoint/windway/voicing geometry | Existing physical headjoint |
| Standard whistle fingering chart | Spec, fixed in v1 |
| Hole positions, diameters, acoustic corrections | Acoustic optimizer |
| Fold path meeting low A length and finger layout | Fold-path solver |
| Exterior body aesthetics and fold topology choices | LLM + human, iterated in build123d-mcp |
| Incremental CAD modeling, measuring, rendering, export | LLM using build123d-mcp |
| Chimney depths, clearances, split shell, STL | Geometry engine through build123d-mcp |
| Final say on tuning and speaking quality | Human with tuner |

## 5. Data Model - `InstrumentSpec`

One JSONC document, with a Python dataclass mirror. All spec lengths are in
millimetres unless the field name says otherwise. Frequencies are in Hz. The
Openwind adapter converts to SI units at the boundary.

The acoustic arc-length origin is the headjoint window datum. The printable body
usually starts downstream at the headjoint socket. Therefore the full acoustic
bore length is:

```
headjoint.window_to_body_socket + body_centerline.arc_length
```

### Example Seed Spec

This example is intentionally a seed, not a final solved instrument. The user
must measure the existing headjoint and the low A reference body length/hole
layout before M1.

```jsonc
{
  "name": "low-d-folded-whistle-v1",
  "target": {
    "instrument": "low_d_tin_whistle",
    "fundamental": "D4",
    "fundamental_hz": 293.66,
    "pitch_standard_hz": 440.0,
    "temperament": "equal",
    "required_notes": [
      "D4", "E4", "F#4", "G4", "A4", "B4", "C#5", "D5",
      "E5", "F#5", "G5", "A5", "B5", "C#6", "D6"
    ]
  },

  "headjoint": {
    "reference": "reference-headjoints/user-low-d-headjoint-v1.step",
    "fixed": true,
    "acoustic_origin": "window_center",
    "measurements_confirmed": true,
    "tube_outer_diameter": 25.0,
    "tube_inner_diameter": 22.0,
    "tube_inner_radius": 11.0,
    "window_to_tube_end": 120.0,
    "sleeve_extension_beyond_tube_end": 56.0,
    "sleeve_end_from_window": 176.0,
    "sleeve_inner_diameter": 25.0,
    "sleeve_outer_diameter": 28.6,
    "sleeve_fit_clearance_radial": 0.1,
    "sleeve_fit_clearance_diametral": 0.2,
    "printed_body_mating_outer_diameter": 24.8,
    "body_connection": "printed_body_outer_diameter_slips_into_25mm_id_sleeve",
    "window_to_body_socket": 120.0,
    "socket_outer_diameter": 25.0,
    "socket_inner_radius": 11.0,
    "socket_depth": null,
    "nominal_insertion_depth": null,
    "effective_length_correction": 0.0,
    "notes": "Measured from the existing physical headjoint. The printed body acoustic start is at the 120 mm tube end; the 56 mm sleeve is mechanical overlap/retention rather than extra acoustic bore. CAD must not modify voicing geometry."
  },

  "bore": {
    "profile": [
      [0.0, 11.0],
      [120.0, 11.0],
      [590.0, 11.0]
    ],
    "cross_section": "circular",
    "prototype_bore_choice": "cylindrical",
    "body_start_arc": 120.0
  },

  "centerline": {
    "description": "Folded path for the printable body only. It starts at body_start_arc.",
    "segments": [
      {"type": "straight", "length": 172.168},
      {"type": "bend", "angle_deg": 180.0, "radius": 40.0},
      {"type": "straight", "length": 172.168}
    ]
  },

  "holes": [
    {
      "id": "h1",
      "label": "top",
      "finger": "L1",
      "position": null,
      "position_seed": 275.0,
      "diameter": null,
      "diameter_seed": 8.0,
      "position_bounds": [210.0, 340.0],
      "diameter_bounds": [5.5, 11.5],
      "chimney_height_estimate": 3.0,
      "chimney_height_measured": null,
      "exit_normal": [0.0, 0.0, 1.0]
    },
    {
      "id": "h2",
      "label": "second",
      "finger": "L2",
      "position": null,
      "position_seed": 325.0,
      "diameter": null,
      "diameter_seed": 8.0,
      "position_bounds": [260.0, 390.0],
      "diameter_bounds": [5.5, 11.5],
      "chimney_height_estimate": 3.0,
      "chimney_height_measured": null,
      "exit_normal": [0.0, 0.0, 1.0]
    },
    {
      "id": "h3",
      "label": "third",
      "finger": "L3",
      "position": null,
      "position_seed": 380.0,
      "diameter": null,
      "diameter_seed": 8.0,
      "position_bounds": [310.0, 445.0],
      "diameter_bounds": [5.5, 11.5],
      "chimney_height_estimate": 3.0,
      "chimney_height_measured": null,
      "exit_normal": [0.0, 0.0, 1.0]
    },
    {
      "id": "h4",
      "label": "fourth",
      "finger": "R1",
      "position": null,
      "position_seed": 425.0,
      "diameter": null,
      "diameter_seed": 8.0,
      "position_bounds": [360.0, 500.0],
      "diameter_bounds": [5.5, 12.0],
      "chimney_height_estimate": 3.0,
      "chimney_height_measured": null,
      "exit_normal": [0.0, 0.0, 1.0]
    },
    {
      "id": "h5",
      "label": "fifth",
      "finger": "R2",
      "position": null,
      "position_seed": 475.0,
      "diameter": null,
      "diameter_seed": 8.5,
      "position_bounds": [405.0, 545.0],
      "diameter_bounds": [5.5, 12.5],
      "chimney_height_estimate": 3.0,
      "chimney_height_measured": null,
      "exit_normal": [0.0, 0.0, 1.0]
    },
    {
      "id": "h6",
      "label": "bottom",
      "finger": "R3",
      "position": null,
      "position_seed": 525.0,
      "diameter": null,
      "diameter_seed": 9.0,
      "position_bounds": [455.0, 575.0],
      "diameter_bounds": [5.5, 13.0],
      "chimney_height_estimate": 3.0,
      "chimney_height_measured": null,
      "exit_normal": [0.0, 0.0, 1.0]
    }
  ],

  "fingering_scope": "standard six-hole tin-whistle diatonic fingerings, sounding D -> D' -> D''",
  "fingerings": [
    {"note": "D4",  "register": 1, "holes": {"h1": "closed", "h2": "closed", "h3": "closed", "h4": "closed", "h5": "closed", "h6": "closed"}},
    {"note": "E4",  "register": 1, "holes": {"h1": "closed", "h2": "closed", "h3": "closed", "h4": "closed", "h5": "closed", "h6": "open"}},
    {"note": "F#4", "register": 1, "holes": {"h1": "closed", "h2": "closed", "h3": "closed", "h4": "closed", "h5": "open",   "h6": "open"}},
    {"note": "G4",  "register": 1, "holes": {"h1": "closed", "h2": "closed", "h3": "closed", "h4": "open",   "h5": "open",   "h6": "open"}},
    {"note": "A4",  "register": 1, "holes": {"h1": "closed", "h2": "closed", "h3": "open",   "h4": "open",   "h5": "open",   "h6": "open"}},
    {"note": "B4",  "register": 1, "holes": {"h1": "closed", "h2": "open",   "h3": "open",   "h4": "open",   "h5": "open",   "h6": "open"}},
    {"note": "C#5", "register": 1, "holes": {"h1": "open",   "h2": "open",   "h3": "open",   "h4": "open",   "h5": "open",   "h6": "open"}},
    {"note": "D5",  "register": 2, "holes": {"h1": "closed", "h2": "closed", "h3": "closed", "h4": "closed", "h5": "closed", "h6": "closed"}},
    {"note": "E5",  "register": 2, "holes": {"h1": "closed", "h2": "closed", "h3": "closed", "h4": "closed", "h5": "closed", "h6": "open"}},
    {"note": "F#5", "register": 2, "holes": {"h1": "closed", "h2": "closed", "h3": "closed", "h4": "closed", "h5": "open",   "h6": "open"}},
    {"note": "G5",  "register": 2, "holes": {"h1": "closed", "h2": "closed", "h3": "closed", "h4": "open",   "h5": "open",   "h6": "open"}},
    {"note": "A5",  "register": 2, "holes": {"h1": "closed", "h2": "closed", "h3": "open",   "h4": "open",   "h5": "open",   "h6": "open"}},
    {"note": "B5",  "register": 2, "holes": {"h1": "closed", "h2": "open",   "h3": "open",   "h4": "open",   "h5": "open",   "h6": "open"}},
    {"note": "C#6", "register": 2, "holes": {"h1": "open",   "h2": "open",   "h3": "open",   "h4": "open",   "h5": "open",   "h6": "open"}},
    {"note": "D6",  "register": 3, "holes": {"h1": "closed", "h2": "closed", "h3": "closed", "h4": "closed", "h5": "closed", "h6": "closed"}}
  ],

  "physical_layout": {
    "reference": "mk-pro-a-web-derived",
    "source_url": "https://mkwhistles.com/products/mk-pro-low-a-whistle",
    "scale_drawing_url": "https://mkwhistles.com/cdn/shop/files/Pro_A_Scaling_Drawing_v2_7549fd57-1002-4cb1-be67-190358c345fa.pdf",
    "source_confidence": "approximate; derived from published scale drawing using its 20 mm calibration square",
    "meaning": "Physical A whistle body length and approximate finger-hole positions only; acoustics remain low D.",
    "coordinate_system": "x from body socket toward foot, y lateral, z outward from nominal front surface",
    "overall_reference_length": 386.0,
    "body_length": 286.4,
    "body_length_tolerance": 10.0,
    "non_acoustic_tail_allowed": false,
    "length_conflict_policy": "relax_visual_length",
    "max_width": null,
    "max_depth": null,
    "hole_centers_from_socket": {
      "h1": 91.3,
      "h2": 116.8,
      "h3": 143.3,
      "h4": 177.3,
      "h5": 197.3,
      "h6": 230.8
    },
    "hole_position_tolerance": 8.0
  },

  "ergonomics": {
    "finger_template_source": "measured_low_a_body_layout",
    "finger_template": {
      "h1": {"target_xyz": null, "surface_normal": [0.0, 0.0, 1.0], "tolerance": 4.0},
      "h2": {"target_xyz": null, "surface_normal": [0.0, 0.0, 1.0], "tolerance": 4.0},
      "h3": {"target_xyz": null, "surface_normal": [0.0, 0.0, 1.0], "tolerance": 4.0},
      "h4": {"target_xyz": null, "surface_normal": [0.0, 0.0, 1.0], "tolerance": 4.0},
      "h5": {"target_xyz": null, "surface_normal": [0.0, 0.0, 1.0], "tolerance": 4.0},
      "h6": {"target_xyz": null, "surface_normal": [0.0, 0.0, 1.0], "tolerance": 4.0}
    }
  },

  "body": {
    "style_prompt": "low A-sized folded low D whistle body",
    "min_wall": 2.0,
    "shell_split_plane": "bore_plane",
    "split_strategy": "design_split_from_start_with_single_canonical_solid",
    "shell_join": "glued_centerline_seam",
    "gasket": false,
    "alignment_features": true
  },

  "fabrication": {
    "process": "FDM",
    "material": "PETG",
    "print_bed": [256.0, 256.0],
    "max_exported_part_footprint": [250.0, 250.0],
    "layer_orientation": "split halves printed from the centreline split layout, subject to CAD printability checks",
    "shell_closure": "glued_down_centreline",
    "gasket_material": null,
    "headjoint_fit": "push_fit_printed_spigot_into_25mm_id_sleeve",
    "fallback_headjoint_seal": "waxed_thread_wrap_if_push_fit_leaks_or_slips",
    "printability": {
      "mcp_tool": "build123d-mcp.analyze_printability",
      "backend": "pzfreo/augura",
      "analysis_type": "BREP-exact FDM printability",
      "run_on": ["left_shell", "right_shell"],
      "build_volume": [256.0, 256.0, null],
      "bed_fit_margin": 3.0,
      "required_checks": [
        "watertight_manifold",
        "overhangs",
        "thin_walls",
        "minimum_vertical_feature",
        "tip_over_risk",
        "brim_or_raft_need",
        "bed_fit"
      ],
      "support_policy": "prefer redesign/orientation changes over supports inside toneholes, bore, socket, or glue lands"
    }
  },

  "calibration": {
    "temperature_c": 20.0,
    "global_length_correction": 0.0,
    "headjoint_effective_length_correction": 0.0,
    "hole_radius_corrections": {},
    "chimney_corrections": {},
    "observed_fingering_residuals_cents": {}
  }
}
```

### Invariants

- `bore.body_start_arc == headjoint.window_to_body_socket`.
- `centerline.arc_length == bore.profile[-1][0] - bore.body_start_arc` within
  0.25 mm. Bend segment arc length is `abs(angle_rad) * radius`.
- A loaded seed spec may have `holes[].position == null` and
  `holes[].diameter == null`; an optimized spec must not.
- Every optimized hole position must be greater than `bore.body_start_arc` and
  less than the final bore arc length.
- `chimney_height_estimate` may be used only before the first CAD build.
  `chimney_height_measured` must be written back from the actual CAD body before
  convergence is accepted.
- `fingerings` is the source of truth for which holes are open or closed for a
  target note. The LLM does not infer fingering masks from note names.
- The headjoint reference geometry is immutable. CAD tools may model sockets,
  adapters, or mating surfaces, but must not modify windway/window/blade
  geometry.
- `physical_layout` is a low A visual/ergonomic reference only. It must never
  change the target pitch from D4.
- A hole target is geometrically reachable only if the distance from the body
  socket to the target point is less than or equal to that hole's available
  acoustic body arc length, within tolerance:
  `distance(socket, target_xyz) <= hole.position - bore.body_start_arc + tolerance`.
  A fold can make an acoustic path longer than the straight-line distance, but it
  cannot make a point farther away than the arc length available to it.
- `physical_layout.non_acoustic_tail_allowed == false` for this draft. If low A
  visual length conflicts with low D acoustic reach, relax the visual length
  rather than adding a non-acoustic tail.

## 6. Acoustic Engine

Openwind runs in the same Python session as the build123d MCP work so the spec,
CAD-derived chimney depths, and acoustic helpers share one namespace. Use the
project `.mcp.json` to launch `build123d-mcp` via `uv tool run`. The effective
command should be equivalent to:

```
uv tool run --python 3.12 \
  --with openwind \
  --with scipy \
  build123d-mcp@latest \
  --allow-imports openwind,scipy,numpy
```

At session start, call the MCP `version` tool and confirm `build123d-mcp` and
`augura` are listed, then call `health_check`. Augura is the printability backend
used by `build123d-mcp.analyze_printability`; do not replace that gate with a
visual inspection or slicer screenshot.

If Openwind's transitive imports are blocked by the sandbox, use an OS-level
sandbox or local trusted environment and relax the Python import allowlist there.
Do not vendor Openwind into build123d-mcp.

### 6.1 LLM build123d-mcp Workflow

The LLM must use `build123d-mcp` as the normal CAD workbench, not only as a final
export tool. The expected loop for geometry work is:

1. Start each design session with `version`, `health_check`, and `session_state`.
   Confirm `build123d-mcp`, `augura`, and required CAD dependencies are available.
2. Build geometry incrementally with `execute`, keeping stable named parameters
   for all important dimensions: bore radius, bend radius, sleeve fit, wall
   thickness, split-plane offsets, glue-land dimensions, and print-bed limits.
3. Register important shapes with `show(..., name)`. At minimum use names for
   `air_solid`, `body_solid`, `left_shell`, `right_shell`, `headjoint_reference`
   or `headjoint_socket`, and each diagnostic clearance body if created.
4. After every meaningful boolean or sweep, call `measure` before relying on a
   render. Numeric checks are the first proof that the model changed correctly.
5. Use `render_view` for human and LLM visual review after `measure` confirms the
   object is plausible. Use clip planes to inspect bore continuity, toneholes,
   glue lands, and socket geometry.
6. Use `clearance`, `cross_sections`, `find_holes`, `resolve`, and
   `verify_spec` for hard geometric checks: adjacent bore passes, chimneys vs.
   bore passes, shell wall thickness, sleeve fit, and print-bed footprint.
7. Use `save_snapshot` before risky edits and `diff_snapshot` or
   `shape_compare` after changes that should affect only a local region.
8. Use `script` or an equivalent reproducible artifact before export. The final
   STL/STEP files must be reproducible from the saved script plus
   `InstrumentSpec`.
9. Use `export` only after hard checks pass: validity, dimensions, clearances,
   chimney write-back, and Augura printability.

The LLM may propose topology and body-shape changes, but it must not accept an
unmeasured visual impression as proof. Any unresolved MCP error, failed boolean,
unexpected topology count, or failed printability report sends the design back to
the relevant solver/CAD step.

Core functions in the foldflute prelude:

- `load_spec(data) -> InstrumentSpec` validates units, invariants, and the
  fingering chart.
- `to_openwind_geometry(spec, fingering) -> OpenwindInput` converts mm to metres,
  applies headjoint effective-length correction, applies bend corrections, and
  uses measured chimney heights when available.
- `analyze(spec) -> {note: {resonance_hz, cents_vs_target, mode_index}}`
  computes input impedance per fingering and reads the impedance **minima**
  relevant to open-open flute/whistle resonators.
- `optimize_holes(spec) -> spec` solves tonehole positions and diameters for the
  required standard fingerings. It uses Openwind inversion if the constraint set
  is expressible there; otherwise it falls back to scipy optimization around
  `analyze`.
- `apply_calibration(spec, observations) -> spec` fits physical corrections:
  global length, headjoint effective length, hole-radius offsets, and chimney
  offsets. Per-fingering cents residuals are stored for diagnostics, not used as
  the primary fitted variables.

### 6.2 Bend Correction Layer

Openwind has no native 3D bend model. The spec-to-openwind translation applies
published effective-length corrections for bent air columns per bend, as a
function of bend angle and radius/bore ratio. Corrections live in one cited module
so they can be refined against calibration data.

Loss and radiation defaults:

- viscothermal losses on;
- unflanged radiation at open toneholes and the foot;
- room temperature from `spec.calibration.temperature_c`;
- no undercutting model in v1.

## 7. Fold-Path Solver

Input:

- optimized acoustic hole positions and diameters from section 6;
- body centerline topology: fold count, segment order, bend angles, bend radii;
- low A physical layout: body length and approximate hole centers from the socket;
- optional width/depth envelope limits;
- per-hole target exterior point, tolerance, and preferred surface normal.

Output:

- centerline parameters such that each acoustic hole position exits the body
  within tolerance of its assigned low A-style finger target, where reachable;
- hole exit axes that are sealable by the assigned finger and do not intersect
  another bore pass;
- an infeasibility report if the requested topology cannot satisfy the constraints.

Method:

- constrained least squares over segment lengths, bend radii, bend orientations,
  and optional body pose;
- hard rejection for wall, clearance, bend-radius, and physical reachability
  violations;
- soft penalties for finger target error and deviation from preferred normals.

If infeasible, the solver reports which hole/finger/layout constraint failed. The
LLM changes topology or relaxes the low A visual layout; it does not add
non-acoustic tail length or hand-edit acoustic hole positions.

## 8. Geometry Engine and Verification

Build steps, all inside the `build123d-mcp` session and following the section
6.1 workflow:

1. Import or parametrically model the existing headjoint interface. Treat the
   voicing geometry as immutable reference.
2. Build one canonical folded body centerline beginning at `bore.body_start_arc`.
3. Sweep the circular bore along the centerline to create one canonical air solid.
4. Choose the centreline split plane, glue lands, alignment features, and print
   orientation before finalizing the exterior body.
5. Add socket/adapter geometry to mate with the existing headjoint.
6. Generate the exterior body to match the low A reference length and approximate
   hole layout where acoustically reachable.
7. Boolean-subtract the air solid.
8. Cut tonehole chimneys from bore wall to body surface at solved acoustic
   positions and solved exit normals.
9. Write back actual `chimney_height_measured` for every hole.
10. Register the two split shell solids as named MCP objects, conventionally
    `left_shell` and `right_shell`.
11. Run `analyze_printability` on each shell object using the Augura backend and
    store both the plain-text summary and JSON findings with the build artifacts.
12. Export the two split shells and a reproducible build script/artifact only
    after the Augura-backed printability gate passes.

Hard constraints:

- Bend radius >= 1.5x local bore diameter; 2x preferred.
- Folded body approximates the low A physical length within tolerance where
  compatible with low D acoustic reach. If length conflicts, relax visual length;
  do not add a non-acoustic tail in this draft.
- Hole exits are within tolerance of the low A reference positions, unless the
  solver has reported the exact low A location as geometrically unreachable and
  the layout has been relaxed.
- Wall between adjacent bore passes >= `body.min_wall`.
- No tonehole chimney intersects another bore pass.
- Glued centreline seam has continuous bonding land and enough access for clean
  assembly.
- Each exported shell fits inside the 256 mm x 256 mm print bed with margin;
  default maximum footprint is 250 mm x 250 mm.
- Headjoint socket fit has explicit tolerance and insertion stop. The baseline is
  push-fit into the 25 mm ID sleeve; waxed thread wrap is allowed as a fallback
  seal/retention method if the printed fit leaks or slips.
- `build123d-mcp.analyze_printability` passes on every shell, using Augura
  (`pzfreo/augura`) as the printability backend.
- No moisture trap: the lowest interior point in playing orientation drains
  toward an opening, or a drain feature is added.
- Textured or sculptural body details must not thin the wall below `min_wall`.

### 8.1 MCP / Augura Printability Protocol

Printability analysis is an MCP-driven hard gate, not a manual review. For every
candidate geometry:

1. Use `show(shell, "left_shell")` and `show(shell, "right_shell")` in
   `build123d-mcp` so each half is addressable by name.
2. Run `measure` on both shells to confirm non-zero volume, expected bounding
   boxes, and PETG material mass estimates.
3. Run `clearance` and local wall checks for bore-to-bore and chimney-to-bore
   spacing before printability analysis.
4. Run `analyze_printability` on `left_shell` and `right_shell`, with build
   volume set to the 256 mm x 256 mm bed and the default 250 mm x 250 mm maximum
   footprint margin from `spec.fabrication`.
5. Treat any Augura FAIL as design rejection. Preferred fixes are fold topology,
   shell split/orientation, glue land geometry, wall thickness, or local body
   shape changes. Do not accept supports inside the bore, toneholes, socket, or
   glue lands as a primary fix.
6. Store the Augura JSON report beside the STL/STEP exports and include the
   pass/fail summary in the build notes.

Required Augura-backed findings for each shell:

- watertight/manifold solid;
- overhangs compatible with the chosen FDM orientation, or explicitly redesigned;
- thin walls no thinner than `body.min_wall`;
- minimum vertical feature compatible with the intended layer height;
- tip-over/brim/raft risk acceptable for PETG on the selected bed;
- bed fit inside 256 mm x 256 mm, with exported footprint no larger than
  250 mm x 250 mm unless intentionally overridden in the spec.

## 9. LLM Orchestration

The LLM drives the loop through the `build123d-mcp` session. Deterministic
helpers enforce the invariants; the LLM proposes topology, low A layout
interpretation, and body aesthetics. Every geometry iteration should leave MCP
evidence: named session objects, measurements, renders or clipped renders,
clearance checks, snapshots/diffs where useful, and printability reports before
export.

1. Human supplies measured headjoint dimensions or a headjoint CAD/scan.
2. Human supplies the low A reference body length and approximate hole/finger
   positions.
3. LLM drafts `InstrumentSpec` with the standard whistle fingering chart.
4. `optimize_holes` solves acoustic positions and diameters.
5. Fold-path solver places those positions into the low A-style physical layout
   where reachable.
6. Geometry build in `build123d-mcp` writes back measured chimney depths.
7. `analyze_printability` runs through build123d-mcp/Augura on both shell
   objects; rejected geometry returns to topology/body/split design.
8. `analyze` reruns with measured chimney depths and bend corrections.
9. Iterate steps 4-8 until cents table, hole dimensions, chimney depths, and
   Augura printability reports are stable.
10. Render for human review, export STLs, print prototype.

The LLM never hand-adjusts hole positions, hole diameters, or bend corrections.
It changes topology, layout choices, body shape, or targets and reruns the
solvers.

## 10. Calibration Protocol

1. Glue the PETG split shells along the centreline seam and let the adhesive cure.
2. Assemble the existing headjoint with the printed prototype body at the nominal
   insertion depth using the sleeve push fit.
3. Leak-test the glued shell and sleeve joint: block all holes and the foot, then
   check for back-pressure. If the sleeve joint leaks or slips, wrap the mating
   spigot with waxed thread and repeat the leak test.
4. Record each required fingering with a tuner: three trials, moderate breath,
   stable headjoint insertion, temperature noted.
5. Fit physical corrections: global body length, headjoint effective length,
   hole-radius offsets, and chimney offsets.
6. Store observations and fitted corrections in `spec.calibration`.
7. Re-run optimizer and geometry, then reprint if corrections exceed tolerances.

Acceptance after one calibration iteration:

- D4 through D6 standard whistle fingerings within +/-10 cents at moderate breath;
- no hole requires uncomfortable half-covering for a required diatonic note;
- headjoint socket remains airtight and repeatable after disassembly/reassembly.

## 11. Milestones

- **M0 - Headjoint capture.** Measure the existing headjoint: socket dimensions,
  outlet bore radius, window-to-body-socket distance, nominal insertion depth,
  and source/provenance. Exit: reference CAD or parametric interface exists and
  the spec can validate against it.
- **M1 - Straight-body low D validation.** Print a straight low D body for the
  existing headjoint with standard whistle fingerings. Exit: toolchain runs
  end-to-end and relative tuning is within +/-15 cents before calibration.
- **M2 - Folded low D in low A physical layout.** Fold the same acoustic design
  so the exterior length and approximate hole positions follow the measured low A
  reference where reachable. Exit: folded version is within +/-10 cents of the
  straight sibling per required fingering after applying the same calibration.
- **M3 - Ergonomic shell.** Refine the exterior body and finger template for the
  target hand/body size. Exit: all required holes are reachable and sealable,
  and build123d-mcp/Augura printability checks pass for both split shells.
- **M4 - LLM end-to-end.** LLM drives measured brief -> solved spec -> validated
  print files. Exit: one printed folded low D whistle accepted by the human
  player.

## 12. Licensing

- build123d-mcp is Apache-2.0.
- Augura (`pzfreo/augura`) is Apache-2.0 and is used through
  `build123d-mcp.analyze_printability` for BREP-exact FDM printability analysis.
- Openwind is GPLv3 as distributed on PyPI at the time of this draft. Keep it as
  a user-installed optional runtime dependency, for example via `uv tool run
  --with openwind`, and do not make it a bundled dependency of build123d-mcp.
- foldflute imports Openwind in its acoustic helpers. If distributed, license
  foldflute GPL-compatibly or split the Openwind-dependent pieces into a
  separately licensed/private component.
- The existing headjoint CAD/scan must have clear provenance. If it is derived
  from a commercial whistle, keep it private unless permission is explicit.

## 13. Risks and Open Questions

| Risk | Mitigation |
|---|---|
| Existing headjoint effective length is unknown | M0/M1 straight-body calibration before any folded design |
| Headjoint speaks poorly with the solved bore diameter | Match body inlet radius to measured headjoint outlet; test straight body first |
| Fold bend corrections exceed model validity | Keep bend radius >= 2x bore diameter where the layout allows; compare M2 to M1 |
| Deep chimneys from sculptural body detune notes | CAD writes measured chimney heights back before acoustic convergence |
| Split-shell leaks alter tuning | Continuous centreline glue land, alignment features, adhesive cure time, and leak test before tuning measurement |
| Sleeve push fit leaks or slips | Target 24.8 mm printed spigot OD for 25 mm ID sleeve; add waxed thread wrap if needed |
| Augura rejects a shell for overhang, thin wall, bed fit, or tip-over risk | Change split orientation, body shape, glue lands, fold topology, or wall thickness; do not rely on supports inside acoustic or mating features |
| Low A hole targets are farther from the socket than low D acoustic arc length permits | Solver reports infeasible hole/finger pairs; relax the visual layout rather than adding non-acoustic tail length |
| PETG print material, temperature, or moisture drift | Record PETG brand, layer orientation, temperature, and headjoint insertion in calibration |

Open questions:

None for draft v0.3.
