from build123d import *

import numpy, scipy, openwind
from openwind import ImpedanceComputation
print("imports OK: numpy, scipy, openwind, ImpedanceComputation")

import numpy as np
from openwind import ImpedanceComputation

fs = np.linspace(200.0, 2000.0, 20)
geom = [[0.0, 0.011], [0.59, 0.011]]

holes_v1 = [['label', 'position', 'radius', 'chimney'],
            ['h1', 0.4, 0.004, 0.003]]
chart = [['label', 'NoteA'], ['h1', 'x']]
try:
    r = ImpedanceComputation(fs, geom, holes_v1, chart, note='NoteA', temperature=20.0)
    print('variant1 OK')
except Exception as e:
    print('variant1 fail:', type(e).__name__, str(e)[:300])

holes_v2 = [['x', 'l', 'r', 'label'],
            [0.4, 0.003, 0.004, 'h1']]
try:
    r = ImpedanceComputation(fs, geom, holes_v2, chart, note='NoteA', temperature=20.0)
    print('variant2 OK')
except Exception as e:
    print('variant2 fail:', type(e).__name__, str(e)[:300])

import numpy as np
from openwind import ImpedanceComputation

fs = np.linspace(200.0, 2000.0, 40)
geom = [[0.0, 0.011], [0.59, 0.011]]
holes = [['label', 'position', 'radius', 'chimney'],
         ['h1', 0.4, 0.004, 0.003]]
chart = [['label', 'NoteA'], ['h1', 'x']]
r = ImpedanceComputation(fs, geom, holes, chart, note='NoteA', temperature=20.0)

names = [n for n in dir(r) if not n.startswith('_')]
print(sorted(names))
print('Zc?', hasattr(r, 'Zc'))
Z = r.impedance
print('impedance shape:', Z.shape, 'dtype:', Z.dtype)
print('|Z/Zc| first 5:', np.abs(Z / r.Zc)[:5])

# antiresonance helper?
if hasattr(r, 'antiresonance_frequencies'):
    try:
        print('antires:', r.antiresonance_frequencies())
    except Exception as e:
        print('antires call fail:', str(e)[:200])

# radiation kwarg probe
try:
    r2 = ImpedanceComputation(fs, geom, holes, chart, note='NoteA', temperature=20.0,
                              radiation_category='unflanged')
    print('radiation_category unflanged OK')
except Exception as e:
    print('radiation kwarg fail:', type(e).__name__, str(e)[:200])

import numpy as np
from openwind import ImpedanceComputation

geom = [[0.0, 0.011], [0.59, 0.011]]
chart = [['label', 'closed', 'open'], ['h1', 'x', 'o']]
fs = np.linspace(250.0, 500.0, 60)

def first_min_hz(r):
    a = np.abs(r.impedance)
    i = int(np.argmin(a))
    return r.frequencies[i]

# variant1: label/position/radius/chimney -> radius=0.004, chimney=0.003
holes_v1 = [['label', 'position', 'radius', 'chimney'], ['h1', 0.4, 0.004, 0.003]]
r1 = ImpedanceComputation(fs, geom, holes_v1, chart, note='open', temperature=20.0)
print('v1 open-hole first |Z| min at', round(first_min_hz(r1), 1), 'Hz')

# variant2: x/l/r/label -> l(chimney)=0.003, r=0.004  (should match v1 if consistent)
holes_v2 = [['x', 'l', 'r', 'label'], [0.4, 0.003, 0.004, 'h1']]
r2 = ImpedanceComputation(fs, geom, holes_v2, chart, note='open', temperature=20.0)
print('v2 open-hole first |Z| min at', round(first_min_hz(r2), 1), 'Hz')

# closed for reference
r1.set_note('closed')
r1.recompute_impedance_at(np.linspace(250.0, 320.0, 40))
print('after set_note+recompute: n freqs =', len(r1.frequencies), 'first min at',
      round(first_min_hz(r1), 1), 'Hz')

# swap radius/chimney in v1 to prove column meaning differs
holes_v1b = [['label', 'position', 'radius', 'chimney'], ['h1', 0.4, 0.003, 0.004]]
r3 = ImpedanceComputation(fs, geom, holes_v1b, chart, note='open', temperature=20.0)
print('v1 swapped r/chimney open min at', round(first_min_hz(r3), 1), 'Hz')

# ================= foldflute prelude 1/3: spec + acoustics =================
# Units: mm and Hz at all public boundaries; SI (m) only inside the openwind
# adapter (SPEC.md section 5/6). All public functions use plain dicts.
import math
import copy
import numpy as np
from openwind import ImpedanceComputation

# ---------------- notes ----------------
NOTE_SEMITONES = {'C': -9, 'C#': -8, 'Db': -8, 'D': -7, 'D#': -6, 'Eb': -6,
                  'E': -5, 'F': -4, 'F#': -3, 'Gb': -3, 'G': -2, 'G#': -1,
                  'Ab': -1, 'A': 0, 'A#': 1, 'Bb': 1, 'B': 2}

def note_hz(name, pitch_standard_hz=440.0):
    letter = name[:-1]
    octave = int(name[-1])
    semis = NOTE_SEMITONES[letter] + 12 * (octave - 4)
    return pitch_standard_hz * (2.0 ** (semis / 12.0))

def cents_of(f, f_ref):
    return 1200.0 * math.log(f / f_ref) / math.log(2.0)

# ---------------- spec helpers ----------------
def segment_arc(seg):
    if seg['type'] == 'straight':
        return seg['length']
    if seg['type'] == 'bend':
        return abs(math.radians(seg['angle_deg'])) * seg['radius']
    raise ValueError('unknown centerline segment type: ' + str(seg['type']))

def centerline_arc_length(centerline):
    return sum(segment_arc(s) for s in centerline['segments'])

def bore_end(spec):
    return spec['bore']['profile'][-1][0]

def bore_radius(spec):
    return spec['bore']['profile'][0][1]

def load_spec(data):
    """Validate an InstrumentSpec dict against SPEC.md section 5 invariants."""
    spec = copy.deepcopy(data)
    errs = []
    hj = spec['headjoint']
    if spec['bore']['body_start_arc'] != hj['window_to_body_socket']:
        errs.append('bore.body_start_arc != headjoint.window_to_body_socket')
    arc = centerline_arc_length(spec['centerline'])
    expect = bore_end(spec) - spec['bore']['body_start_arc']
    if abs(arc - expect) > 0.25:
        errs.append('centerline arc %.3f != bore body length %.3f (tol 0.25)'
                    % (arc, expect))
    xs = [p[0] for p in spec['bore']['profile']]
    if sorted(xs) != xs:
        errs.append('bore profile x not monotone')
    ids = [h['id'] for h in spec['holes']]
    if len(set(ids)) != len(ids):
        errs.append('duplicate hole ids')
    start = spec['bore']['body_start_arc']
    end = bore_end(spec)
    for h in spec['holes']:
        pos = h['position'] if h['position'] is not None else h['position_seed']
        lo, hi = h['position_bounds']
        if not (start < lo < hi < end):
            errs.append('%s position_bounds outside body bore' % h['id'])
        if h['position'] is not None and not (start < h['position'] < end):
            errs.append('%s optimized position outside bore' % h['id'])
    chart_notes = [f['note'] for f in spec['fingerings']]
    for n in spec['target']['required_notes']:
        if n not in chart_notes:
            errs.append('required note %s missing from fingerings' % n)
    for f in spec['fingerings']:
        if sorted(f['holes'].keys()) != sorted(ids):
            errs.append('fingering %s does not cover all holes' % f['note'])
        for st in f['holes'].values():
            if st not in ('open', 'closed'):
                errs.append('fingering %s has invalid state %s' % (f['note'], st))
    if errs:
        raise ValueError('spec invalid:\n  ' + '\n  '.join(errs))
    return spec

def effective_holes(spec, positions=None, diameters=None, chimneys=None):
    """Working values per hole: solved fields, else seeds; optional overrides."""
    out = []
    for h in spec['holes']:
        hid = h['id']
        pos = h['position'] if h['position'] is not None else h['position_seed']
        dia = h['diameter'] if h['diameter'] is not None else h['diameter_seed']
        chy = (h['chimney_height_measured']
               if h['chimney_height_measured'] is not None
               else h['chimney_height_estimate'])
        if positions and hid in positions: pos = positions[hid]
        if diameters and hid in diameters: dia = diameters[hid]
        if chimneys and hid in chimneys: chy = chimneys[hid]
        out.append({'id': hid, 'position': pos, 'diameter': dia, 'chimney': chy})
    return out

def set_bore_end(spec, new_end):
    """Move the bore end, keeping the centerline invariant by adjusting the
    final straight segment."""
    old = bore_end(spec)
    delta = new_end - old
    spec['bore']['profile'][-1][0] = new_end
    last = spec['centerline']['segments'][-1]
    if last['type'] != 'straight':
        raise ValueError('final centerline segment must be straight to retune length')
    if last['length'] + delta <= 0:
        raise ValueError('bore end change exceeds final straight segment')
    last['length'] += delta

# ---------------- bend corrections (SPEC.md 6.2) ----------------
# A toroidal bend of centerline radius R and bore radius a behaves, in the
# long-wavelength limit, as a straight duct of the same volume but reduced
# inertance: with B = a/R and a 1/r potential-flow velocity profile across the
# section, M_bend/M_straight = (1 + sqrt(1 - B^2)) / 2  ~=  1 - B^2/4.
# Reduced inertance with unchanged compliance means the bend is acoustically
# SHORTER than its arc length; pitch rises. The correction magnitude depends
# on the bend's position in the standing wave (inertance- vs compliance-
# dominated); kappa=1 applies the full (velocity-antinode) correction and is
# the calibration-refinable knob.
# Refs: C.J. Nederveen, JASA 104(3):1616-1626 (1998); S. Felix, J.-P. Dalmont,
# C.J. Nederveen, JASA 131:4164-4172 (2012).
BEND_KAPPA_DEFAULT = 1.0

def bend_shrink_factor(bore_radius_mm, bend_radius_mm, kappa=BEND_KAPPA_DEFAULT):
    B = bore_radius_mm / bend_radius_mm
    if B >= 1.0:
        raise ValueError('bend radius must exceed bore radius')
    ratio = (1.0 + math.sqrt(1.0 - B * B)) / 2.0
    return kappa * (1.0 - ratio)   # fraction of bend arc removed

def centerline_bend_spans(spec):
    """Bends as [{'s0': arc-from-window at bend start, 'arc', 'radius'}]."""
    s = spec['bore']['body_start_arc']
    out = []
    for seg in spec['centerline']['segments']:
        a = segment_arc(seg)
        if seg['type'] == 'bend':
            out.append({'s0': s, 'arc': a, 'radius': seg['radius']})
        s += a
    return out

def arc_effective_map(spec, kappa=BEND_KAPPA_DEFAULT):
    """Map arc position (mm from window) -> bend-corrected effective position."""
    bends = centerline_bend_spans(spec)
    a_bore = bore_radius(spec)
    def f(s):
        s_eff = s
        for b in bends:
            shrink = bend_shrink_factor(a_bore, b['radius'], kappa)
            covered = min(max(s - b['s0'], 0.0), b['arc'])
            s_eff -= shrink * covered
        return s_eff
    return f

print('prelude 1/3 part A loaded: notes, spec, bends')

# ================= foldflute prelude 1/3 part B: openwind adapter + analyze =================
MM = 0.001

def to_openwind_geometry(spec, positions=None, diameters=None, chimneys=None,
                         L_end=None, apply_bends=True, kappa=BEND_KAPPA_DEFAULT,
                         extra_holes=None):
    """InstrumentSpec (mm) -> (main_bore_m, holes_table_m, chart, temperature_c).
    Applies headjoint effective-length correction, calibration corrections,
    and toroidal bend corrections. extra_holes: list of
    {'id','position','diameter','chimney','state'} always-open/closed extras
    (e.g. a drain hole experiment)."""
    cal = spec['calibration']
    lead = spec['headjoint']['effective_length_correction']
    end = L_end if L_end is not None else bore_end(spec)
    smap = arc_effective_map(spec, kappa) if apply_bends else (lambda s: s)
    end_eff = smap(end) + cal['global_length_correction'] + lead
    r_m = bore_radius(spec) * MM
    main_bore = [[0.0, r_m], [end_eff * MM, r_m]]

    rows = [['label', 'position', 'radius', 'chimney']]
    hole_states = {}   # id -> None (use chart) or fixed state
    for h in effective_holes(spec, positions, diameters, chimneys):
        hid = h['id']
        rad = h['diameter'] / 2.0 + cal['hole_radius_corrections'].get(hid, 0.0)
        chy = h['chimney'] + cal['chimney_corrections'].get(hid, 0.0)
        pos = smap(h['position']) + lead
        rows.append([hid, pos * MM, rad * MM, chy * MM])
        hole_states[hid] = None
    for h in (extra_holes or []):
        pos = smap(h['position']) + lead
        rows.append([h['id'], pos * MM, h['diameter'] / 2.0 * MM, h['chimney'] * MM])
        hole_states[h['id']] = h['state']

    notes = [f['note'] for f in spec['fingerings']]
    chart = [['label'] + notes]
    state_char = {'closed': 'x', 'open': 'o'}
    for f in spec['fingerings']:
        pass
    for hid in hole_states:
        if hole_states[hid] is None:
            row = [hid]
            for f in spec['fingerings']:
                row.append(state_char[f['holes'][hid]])
        else:
            row = [hid] + [state_char[hole_states[hid]]] * len(notes)
        chart.append(row)
    return main_bore, rows, chart, cal['temperature_c']

def _refine_min(freqs, logmag, i):
    """Parabolic vertex through 3 points around index i in (log f, log|Z|)."""
    l = np.log(freqs[i - 1:i + 2])
    y = logmag[i - 1:i + 2]
    d1 = (y[1] - y[0]) / (l[1] - l[0])
    d2 = (y[2] - y[1]) / (l[2] - l[1])
    curv = (d2 - d1) / ((l[2] - l[0]) / 2.0)
    if curv <= 0:
        return freqs[i]
    slope_mid = (y[2] - y[0]) / (l[2] - l[0])
    l_star = (l[0] + l[2]) / 2.0 - slope_mid / curv
    l_star = min(max(l_star, l[0]), l[2])
    return math.exp(l_star)

def _find_impedance_min(comp, target_hz, span_cents=170.0, n_grid=29):
    """Frequency of the |Z| minimum nearest target_hz. Widens the window if the
    minimum sits on a boundary; parabolic refinement + fine pass."""
    lo, hi = -span_cents, span_cents
    for _ in range(5):
        fs = target_hz * 2.0 ** (np.linspace(lo, hi, n_grid) / 1200.0)
        comp.recompute_impedance_at(fs)
        a = np.log(np.abs(comp.impedance))
        i = int(np.argmin(a))
        if i == 0:
            lo -= span_cents * 0.9
        elif i == n_grid - 1:
            hi += span_cents * 0.9
        else:
            f1 = _refine_min(fs, a, i)
            fs2 = f1 * 2.0 ** (np.linspace(-8.0, 8.0, 9) / 1200.0)
            comp.recompute_impedance_at(fs2)
            a2 = np.log(np.abs(comp.impedance))
            j = int(np.argmin(a2))
            if j == 0 or j == 8:
                return f1
            return _refine_min(fs2, a2, j)
    raise RuntimeError('no impedance minimum near %.1f Hz' % target_hz)

def build_computation(spec, positions=None, diameters=None, chimneys=None,
                      L_end=None, apply_bends=True, kappa=BEND_KAPPA_DEFAULT,
                      extra_holes=None):
    main_bore, holes, chart, temp = to_openwind_geometry(
        spec, positions, diameters, chimneys, L_end, apply_bends, kappa, extra_holes)
    notes = [f['note'] for f in spec['fingerings']]
    fs0 = np.array([200.0, 300.0, 400.0])
    return ImpedanceComputation(fs0, main_bore, holes, chart, note=notes[0],
                                temperature=temp, losses=True,
                                radiation_category='unflanged')

def analyze(spec, notes=None, positions=None, diameters=None, chimneys=None,
            L_end=None, apply_bends=True, kappa=BEND_KAPPA_DEFAULT,
            extra_holes=None, comp=None):
    """Passive resonance prediction per fingering: impedance minima (SPEC.md 6).
    Returns {note: {'target_hz', 'f_res', 'cents'}}."""
    std = spec['target']['pitch_standard_hz']
    if comp is None:
        comp = build_computation(spec, positions, diameters, chimneys, L_end,
                                 apply_bends, kappa, extra_holes)
    want = notes if notes is not None else [f['note'] for f in spec['fingerings']]
    out = {}
    for n in want:
        t = note_hz(n, std)
        comp.set_note(n)
        f_res = _find_impedance_min(comp, t)
        out[n] = {'target_hz': t, 'f_res': f_res, 'cents': cents_of(f_res, t)}
    return out

def cents_table(res, order=None):
    keys = order if order else list(res.keys())
    lines = ['note    target(Hz)  res(Hz)   cents']
    for n in keys:
        r = res[n]
        lines.append('%-6s  %9.2f  %8.2f  %+6.1f'
                     % (n, r['target_hz'], r['f_res'], r['cents']))
    return '\n'.join(lines)

print('prelude 1/3 part B loaded: adapter, resonance finder, analyze')

import json
SPEC_JSON = r'''
{
  "name": "low-d-folded-whistle-v1",
  "target": {
    "instrument": "low_d_tin_whistle",
    "fundamental": "D4",
    "fundamental_hz": 293.66,
    "pitch_standard_hz": 440.0,
    "temperament": "equal",
    "required_notes": ["D4","E4","F#4","G4","A4","B4","C#5","D5","E5","F#5","G5","A5","B5","C#6","D6"]
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
    "notes": "Measured from the existing physical headjoint."
  },
  "bore": {
    "profile": [[0.0, 11.0],[120.0, 11.0],[590.0, 11.0]],
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
    {"id":"h1","label":"top","finger":"L1","position":null,"position_seed":275.0,"diameter":null,"diameter_seed":8.0,"position_bounds":[210.0,340.0],"diameter_bounds":[5.5,11.5],"chimney_height_estimate":3.0,"chimney_height_measured":null,"exit_normal":[0.0,0.0,1.0]},
    {"id":"h2","label":"second","finger":"L2","position":null,"position_seed":325.0,"diameter":null,"diameter_seed":8.0,"position_bounds":[260.0,390.0],"diameter_bounds":[5.5,11.5],"chimney_height_estimate":3.0,"chimney_height_measured":null,"exit_normal":[0.0,0.0,1.0]},
    {"id":"h3","label":"third","finger":"L3","position":null,"position_seed":380.0,"diameter":null,"diameter_seed":8.0,"position_bounds":[310.0,445.0],"diameter_bounds":[5.5,11.5],"chimney_height_estimate":3.0,"chimney_height_measured":null,"exit_normal":[0.0,0.0,1.0]},
    {"id":"h4","label":"fourth","finger":"R1","position":null,"position_seed":425.0,"diameter":null,"diameter_seed":8.0,"position_bounds":[360.0,500.0],"diameter_bounds":[5.5,12.0],"chimney_height_estimate":3.0,"chimney_height_measured":null,"exit_normal":[0.0,0.0,1.0]},
    {"id":"h5","label":"fifth","finger":"R2","position":null,"position_seed":475.0,"diameter":null,"diameter_seed":8.5,"position_bounds":[405.0,545.0],"diameter_bounds":[5.5,12.5],"chimney_height_estimate":3.0,"chimney_height_measured":null,"exit_normal":[0.0,0.0,1.0]},
    {"id":"h6","label":"bottom","finger":"R3","position":null,"position_seed":525.0,"diameter":null,"diameter_seed":9.0,"position_bounds":[455.0,575.0],"diameter_bounds":[5.5,13.0],"chimney_height_estimate":3.0,"chimney_height_measured":null,"exit_normal":[0.0,0.0,1.0]}
  ],
  "fingering_scope": "standard six-hole tin-whistle diatonic fingerings, sounding D -> D' -> D''",
  "fingerings": [
    {"note":"D4","register":1,"holes":{"h1":"closed","h2":"closed","h3":"closed","h4":"closed","h5":"closed","h6":"closed"}},
    {"note":"E4","register":1,"holes":{"h1":"closed","h2":"closed","h3":"closed","h4":"closed","h5":"closed","h6":"open"}},
    {"note":"F#4","register":1,"holes":{"h1":"closed","h2":"closed","h3":"closed","h4":"closed","h5":"open","h6":"open"}},
    {"note":"G4","register":1,"holes":{"h1":"closed","h2":"closed","h3":"closed","h4":"open","h5":"open","h6":"open"}},
    {"note":"A4","register":1,"holes":{"h1":"closed","h2":"closed","h3":"open","h4":"open","h5":"open","h6":"open"}},
    {"note":"B4","register":1,"holes":{"h1":"closed","h2":"open","h3":"open","h4":"open","h5":"open","h6":"open"}},
    {"note":"C#5","register":1,"holes":{"h1":"open","h2":"open","h3":"open","h4":"open","h5":"open","h6":"open"}},
    {"note":"D5","register":2,"holes":{"h1":"closed","h2":"closed","h3":"closed","h4":"closed","h5":"closed","h6":"closed"}},
    {"note":"E5","register":2,"holes":{"h1":"closed","h2":"closed","h3":"closed","h4":"closed","h5":"closed","h6":"open"}},
    {"note":"F#5","register":2,"holes":{"h1":"closed","h2":"closed","h3":"closed","h4":"closed","h5":"open","h6":"open"}},
    {"note":"G5","register":2,"holes":{"h1":"closed","h2":"closed","h3":"closed","h4":"open","h5":"open","h6":"open"}},
    {"note":"A5","register":2,"holes":{"h1":"closed","h2":"closed","h3":"open","h4":"open","h5":"open","h6":"open"}},
    {"note":"B5","register":2,"holes":{"h1":"closed","h2":"open","h3":"open","h4":"open","h5":"open","h6":"open"}},
    {"note":"C#6","register":2,"holes":{"h1":"open","h2":"open","h3":"open","h4":"open","h5":"open","h6":"open"}},
    {"note":"D6","register":3,"holes":{"h1":"closed","h2":"closed","h3":"closed","h4":"closed","h5":"closed","h6":"closed"}}
  ],
  "physical_layout": {
    "reference": "mk-pro-a-web-derived",
    "body_length": 286.4,
    "body_length_tolerance": 10.0,
    "overall_reference_length": 386.0,
    "non_acoustic_tail_allowed": false,
    "length_conflict_policy": "relax_visual_length",
    "max_width": null,
    "max_depth": null,
    "hole_centers_from_socket": {"h1":91.3,"h2":116.8,"h3":143.3,"h4":177.3,"h5":197.3,"h6":230.8},
    "hole_position_tolerance": 8.0
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
    "headjoint_fit": "push_fit_printed_spigot_into_25mm_id_sleeve"
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
'''
spec_seed = load_spec(json.loads(SPEC_JSON))
print('seed spec valid:', spec_seed['name'])

# Straight-body variant for M1: same bore, centerline = one straight segment.
spec_straight = copy.deepcopy(spec_seed)
spec_straight['centerline'] = {
    'description': 'M1 straight body',
    'segments': [{'type': 'straight', 'length': bore_end(spec_seed) - 120.0}]}
spec_straight = load_spec(spec_straight)
print('straight variant valid; body arc =',
      centerline_arc_length(spec_straight['centerline']))

def ow_label(note):
    return note.replace('#', 's')   # '#' is a comment char in openwind charts

def to_openwind_geometry(spec, positions=None, diameters=None, chimneys=None,
                         L_end=None, apply_bends=True, kappa=BEND_KAPPA_DEFAULT,
                         extra_holes=None):
    cal = spec['calibration']
    lead = spec['headjoint']['effective_length_correction']
    end = L_end if L_end is not None else bore_end(spec)
    smap = arc_effective_map(spec, kappa) if apply_bends else (lambda s: s)
    end_eff = smap(end) + cal['global_length_correction'] + lead
    r_m = bore_radius(spec) * MM
    main_bore = [[0.0, r_m], [end_eff * MM, r_m]]

    rows = [['label', 'position', 'radius', 'chimney']]
    hole_states = {}
    for h in effective_holes(spec, positions, diameters, chimneys):
        hid = h['id']
        rad = h['diameter'] / 2.0 + cal['hole_radius_corrections'].get(hid, 0.0)
        chy = h['chimney'] + cal['chimney_corrections'].get(hid, 0.0)
        pos = smap(h['position']) + lead
        rows.append([hid, pos * MM, rad * MM, chy * MM])
        hole_states[hid] = None
    for h in (extra_holes or []):
        pos = smap(h['position']) + lead
        rows.append([h['id'], pos * MM, h['diameter'] / 2.0 * MM, h['chimney'] * MM])
        hole_states[h['id']] = h['state']

    notes = [ow_label(f['note']) for f in spec['fingerings']]
    chart = [['label'] + notes]
    state_char = {'closed': 'x', 'open': 'o'}
    for hid in hole_states:
        if hole_states[hid] is None:
            row = [hid]
            for f in spec['fingerings']:
                row.append(state_char[f['holes'][hid]])
        else:
            row = [hid] + [state_char[hole_states[hid]]] * len(notes)
        chart.append(row)
    return main_bore, rows, chart, cal['temperature_c']

def build_computation(spec, positions=None, diameters=None, chimneys=None,
                      L_end=None, apply_bends=True, kappa=BEND_KAPPA_DEFAULT,
                      extra_holes=None):
    main_bore, holes, chart, temp = to_openwind_geometry(
        spec, positions, diameters, chimneys, L_end, apply_bends, kappa, extra_holes)
    first = chart[0][1]
    fs0 = np.array([200.0, 300.0, 400.0])
    return ImpedanceComputation(fs0, main_bore, holes, chart, note=first,
                                temperature=temp, losses=True,
                                radiation_category='unflanged')

def analyze(spec, notes=None, positions=None, diameters=None, chimneys=None,
            L_end=None, apply_bends=True, kappa=BEND_KAPPA_DEFAULT,
            extra_holes=None, comp=None):
    std = spec['target']['pitch_standard_hz']
    if comp is None:
        comp = build_computation(spec, positions, diameters, chimneys, L_end,
                                 apply_bends, kappa, extra_holes)
    want = notes if notes is not None else [f['note'] for f in spec['fingerings']]
    out = {}
    for n in want:
        t = note_hz(n, std)
        comp.set_note(ow_label(n))
        f_res = _find_impedance_min(comp, t)
        out[n] = {'target_hz': t, 'f_res': f_res, 'cents': cents_of(f_res, t)}
    return out

res = analyze(spec_straight, apply_bends=False)
order = [f['note'] for f in spec_straight['fingerings']]
print(cents_table(res, order))

# ================= foldflute prelude 2/3: hole optimizer =================
# Gauss-Seidel over (bore end -> D4) then (h6..h1 -> E4..C#5), damped secant
# per scalar. Register 2 and D6 are reported as predictions, not targets
# (register balance is set by diameter choices, which are inputs here).
R1_TARGET_NOTES = [('h6', 'E4'), ('h5', 'F#4'), ('h4', 'G4'),
                   ('h3', 'A4'), ('h2', 'B4'), ('h1', 'C#5')]
L_END_BOUNDS = (545.0, 608.0)

def eval_note_cents(spec, note, positions, L_end, diameters=None,
                    chimneys=None, apply_bends=True, kappa=BEND_KAPPA_DEFAULT):
    r = analyze(spec, notes=[note], positions=positions, diameters=diameters,
                chimneys=chimneys, L_end=L_end, apply_bends=apply_bends,
                kappa=kappa)
    return r[note]['cents']

def _secant_solve(fn, x0, bounds, slope0=-2.5, tol=0.5, max_iter=7,
                  max_step=20.0):
    lo, hi = bounds
    x = min(max(x0, lo), hi)
    c = fn(x)
    slope = slope0
    for _ in range(max_iter):
        if abs(c) <= tol:
            break
        step = -c / slope
        step = min(max(step, -max_step), max_step)
        x_new = min(max(x + step, lo), hi)
        if abs(x_new - x) < 5e-4:
            break
        c_new = fn(x_new)
        if abs(x_new - x) > 1e-3:
            s = (c_new - c) / (x_new - x)
            if s < -0.05:
                slope = s
        x, c = x_new, c_new
    return x, c

def optimize_sweep(spec, positions, L_end, diameters=None, chimneys=None,
                   apply_bends=True, kappa=BEND_KAPPA_DEFAULT, tol=0.5):
    """One Gauss-Seidel sweep; mutates positions, returns (positions, L_end)."""
    fnL = lambda L: eval_note_cents(spec, 'D4', positions, L, diameters,
                                    chimneys, apply_bends, kappa)
    L_end, _ = _secant_solve(fnL, L_end, L_END_BOUNDS, slope0=-3.0, tol=tol)
    hole_bounds = {h['id']: tuple(h['position_bounds']) for h in spec['holes']}
    for hid, note in R1_TARGET_NOTES:
        fn = (lambda p, hid=hid, note=note:
              eval_note_cents(spec, note,
                              {**positions, hid: p}, L_end, diameters,
                              chimneys, apply_bends, kappa))
        p_new, _ = _secant_solve(fn, positions[hid], hole_bounds[hid],
                                 slope0=-2.5, tol=tol)
        positions[hid] = p_new
    return positions, L_end

def r1_residuals(spec, positions, L_end, diameters=None, chimneys=None,
                 apply_bends=True, kappa=BEND_KAPPA_DEFAULT):
    notes = ['D4'] + [n for _, n in R1_TARGET_NOTES]
    r = analyze(spec, notes=notes, positions=positions, diameters=diameters,
                chimneys=chimneys, L_end=L_end, apply_bends=apply_bends,
                kappa=kappa)
    return {n: r[n]['cents'] for n in notes}

def writeback_solution(spec, positions, L_end, diameters=None):
    """Solved values -> spec copy (holes[].position/diameter, bore end)."""
    out = copy.deepcopy(spec)
    for h in out['holes']:
        h['position'] = positions[h['id']]
        if diameters and h['id'] in diameters:
            h['diameter'] = diameters[h['id']]
        elif h['diameter'] is None:
            h['diameter'] = h['diameter_seed']
    set_bore_end(out, L_end)
    return load_spec(out)

# ---- initial state for the M1 straight solve (seed diameters) ----
m1_positions = {h['id']: h['position_seed'] for h in spec_straight['holes']}
m1_L = bore_end(spec_straight)
print('optimizer loaded; start positions', m1_positions, 'L_end', m1_L)

m1_positions, m1_L = optimize_sweep(spec_straight, m1_positions, m1_L,
                                    apply_bends=False)
resid = r1_residuals(spec_straight, m1_positions, m1_L, apply_bends=False)
print('after sweep 1: L_end = %.2f' % m1_L)
print({k: round(v, 1) for k, v in m1_positions.items()})
print('R1 residuals (cents):', {k: round(v, 1) for k, v in resid.items()})

def load_spec(data):
    spec = copy.deepcopy(data)
    errs = []
    hj = spec['headjoint']
    if spec['bore']['body_start_arc'] != hj['window_to_body_socket']:
        errs.append('bore.body_start_arc != headjoint.window_to_body_socket')
    arc = centerline_arc_length(spec['centerline'])
    expect = bore_end(spec) - spec['bore']['body_start_arc']
    if abs(arc - expect) > 0.25:
        errs.append('centerline arc %.3f != bore body length %.3f (tol 0.25)'
                    % (arc, expect))
    xs = [p[0] for p in spec['bore']['profile']]
    if sorted(xs) != xs:
        errs.append('bore profile x not monotone')
    ids = [h['id'] for h in spec['holes']]
    if len(set(ids)) != len(ids):
        errs.append('duplicate hole ids')
    start = spec['bore']['body_start_arc']
    end = bore_end(spec)
    for h in spec['holes']:
        lo, hi = h['position_bounds']
        if not (start < lo < hi):
            errs.append('%s position_bounds malformed' % h['id'])
        if h['position'] is not None and not (start < h['position'] < end):
            errs.append('%s optimized position outside bore' % h['id'])
    chart_notes = [f['note'] for f in spec['fingerings']]
    for n in spec['target']['required_notes']:
        if n not in chart_notes:
            errs.append('required note %s missing from fingerings' % n)
    for f in spec['fingerings']:
        if sorted(f['holes'].keys()) != sorted(ids):
            errs.append('fingering %s does not cover all holes' % f['note'])
        for st in f['holes'].values():
            if st not in ('open', 'closed'):
                errs.append('fingering %s has invalid state %s' % (f['note'], st))
    if errs:
        raise ValueError('spec invalid:\n  ' + '\n  '.join(errs))
    return spec

HOLE_END_MARGIN = 8.0   # mm: keep the last hole clear of the open foot

def optimize_sweep(spec, positions, L_end, diameters=None, chimneys=None,
                   apply_bends=True, kappa=BEND_KAPPA_DEFAULT, tol=0.5):
    fnL = lambda L: eval_note_cents(spec, 'D4', positions, L, diameters,
                                    chimneys, apply_bends, kappa)
    L_end, _ = _secant_solve(fnL, L_end, L_END_BOUNDS, slope0=-3.0, tol=tol)
    hole_bounds = {h['id']: (h['position_bounds'][0],
                             min(h['position_bounds'][1],
                                 L_end - HOLE_END_MARGIN))
                   for h in spec['holes']}
    for hid, note in R1_TARGET_NOTES:
        fn = (lambda p, hid=hid, note=note:
              eval_note_cents(spec, note,
                              {**positions, hid: p}, L_end, diameters,
                              chimneys, apply_bends, kappa))
        p_new, _ = _secant_solve(fn, positions[hid], hole_bounds[hid],
                                 slope0=-2.5, tol=tol)
        positions[hid] = p_new
    return positions, L_end

spec_m1 = writeback_solution(spec_straight, m1_positions, m1_L)
res_m1 = analyze(spec_m1, apply_bends=False)
order = [f['note'] for f in spec_m1['fingerings']]
print(cents_table(res_m1, order))
print()
print('solved holes (pos mm from window / dia mm):')
for h in spec_m1['holes']:
    print('  %s  %7.2f  %4.1f' % (h['id'], h['position'], h['diameter']))
print('solved bore end: %.2f mm' % bore_end(spec_m1))

def run_scenario(name, diameters, spec_base, start_positions, start_L,
                 sweeps=2, apply_bends=False, kappa=BEND_KAPPA_DEFAULT):
    pos = dict(start_positions)
    L = start_L
    for i in range(sweeps):
        pos, L = optimize_sweep(spec_base, pos, L, diameters=diameters,
                                apply_bends=apply_bends, kappa=kappa,
                                tol=0.4)
    resid = r1_residuals(spec_base, pos, L, diameters=diameters,
                         apply_bends=apply_bends, kappa=kappa)
    worst_r1 = max(abs(v) for v in resid.values())
    full = analyze(spec_base, positions=pos, diameters=diameters, L_end=L,
                   apply_bends=apply_bends, kappa=kappa)
    return {'name': name, 'diameters': diameters, 'positions': pos,
            'L_end': L, 'worst_r1': worst_r1, 'full': full}

def scenario_summary(sc, spec_base):
    r2 = ['D5', 'E5', 'F#5', 'G5', 'A5', 'B5', 'C#6', 'D6']
    line1 = '%s: L_end=%.1f worst_R1=%.1fc' % (sc['name'], sc['L_end'],
                                               sc['worst_r1'])
    line2 = '  R2 cents: ' + '  '.join('%s %+5.1f' % (n, sc['full'][n]['cents'])
                                       for n in r2)
    span = sc['positions']['h6'] - sc['positions']['h1']
    tgt = spec_base['physical_layout']['hole_centers_from_socket']
    datum = (spec_base['headjoint']['sleeve_end_from_window'])
    errs = {h: (sc['positions'][h] - datum) - tgt[h] for h in tgt}
    line3 = ('  hole span h1-h6 %.1f mm; layout err vs low A (mm): ' % span
             + '  '.join('%s %+.0f' % (h, errs[h]) for h in
                         ['h1', 'h2', 'h3', 'h4', 'h5', 'h6']))
    return '\n'.join([line1, line2, line3])

S1 = {'h1': 11.0, 'h2': 10.5, 'h3': 10.0, 'h4': 10.0, 'h5': 10.5, 'h6': 11.0}
sc1 = run_scenario('S1 big-holes', S1, spec_straight, m1_positions, m1_L)
print(scenario_summary(sc1, spec_straight))

S2 = {'h1': 11.5, 'h2': 11.0, 'h3': 10.5, 'h4': 8.0, 'h5': 7.0, 'h6': 6.5}
sc2 = run_scenario('S2 compress', S2, spec_straight, sc1['positions'],
                   sc1['L_end'])
print(scenario_summary(sc2, spec_straight))
print()
S3 = {'h1': 11.5, 'h2': 11.0, 'h3': 10.5, 'h4': 9.5, 'h5': 9.0, 'h6': 8.5}
sc3 = run_scenario('S3 compromise', S3, spec_straight, sc1['positions'],
                   sc1['L_end'])
print(scenario_summary(sc3, spec_straight))

# ================= foldflute prelude 3/3: fold-path solver =================
# Topology search conclusion (documented in the report): with a 22 mm cylindrical
# bore, a 180-degree bend needs pi*R >= 104 mm of arc (R >= 1.5*bore_d), which fits
# neither between adjacent toneholes (<57 mm gaps), nor above h1 (~80 mm of
# straight available, a fold-back needs >= 207 mm), nor wholly below h6
# (~90 mm tail). The one legal fold is a foot U-bend with h6 a bounded angle
# into the bend. Fold plane is XY (passes side by side laterally); all hole
# exits stay +z as required by the spec exit normals.
BEND_RADIUS_MIN_FACTOR = 1.5
FOOT_TAIL_MM = 6.0          # bore continuing past the bend to the foot opening
H6_MAX_BEND_ANGLE_DEG = 35.0
H5_BEND_CLEARANCE_MM = 12.0
BODY_OD = 28.6              # match sleeve outer diameter
FOOT_CAP_CLEAR = 2.0

def solve_foot_fold(spec, positions, L_end):
    d_bore = 2.0 * bore_radius(spec)
    r_min = BEND_RADIUS_MIN_FACTOR * d_bore
    phi_max = math.radians(H6_MAX_BEND_ANGLE_DEG)
    s5, s6 = positions['h5'], positions['h6']
    body0 = spec['bore']['body_start_arc']
    report = {'feasible': True, 'violations': []}

    R_needed = (L_end - FOOT_TAIL_MM - s6) / (math.pi - phi_max)
    R = max(r_min, R_needed)
    if R_needed < r_min:
        # bend radius clamped up; h6 goes deeper into the bend than phi_max
        report['violations'].append(
            'h6 bend angle exceeds %.0f deg at minimum legal bend radius'
            % H6_MAX_BEND_ANGLE_DEG)
        report['feasible'] = False
    s_bend = L_end - FOOT_TAIL_MM - math.pi * R
    if s_bend < s5 + H5_BEND_CLEARANCE_MM:
        report['violations'].append('bend start intrudes on h5')
        report['feasible'] = False
    phi6 = max(0.0, (s6 - s_bend) / R)

    # exterior positions, body coords: x from body start (arc %.0f), y lateral
    holes_xy = {}
    for hid in positions:
        s = positions[hid]
        if s <= s_bend:
            holes_xy[hid] = (s - body0, 0.0)
        else:
            ph = (s - s_bend) / R
            holes_xy[hid] = ((s_bend - body0) + R * math.sin(ph),
                             R * (1.0 - math.cos(ph)))
    apex_x = (s_bend - body0) + R
    body_len = apex_x + BODY_OD / 2.0 + FOOT_CAP_CLEAR
    foot_x = (s_bend - body0) - FOOT_TAIL_MM   # on the return pass, y = 2R

    # low A layout comparison (SPEC.md physical_layout); datum assumption:
    # x_visible measured from the sleeve end (body emerges), i.e. body x - 56.
    lay = spec['physical_layout']
    datum = (spec['headjoint']['sleeve_end_from_window']
             - spec['headjoint']['window_to_body_socket'])
    tol = lay['hole_position_tolerance']
    hole_err = {}
    for hid, tgt in lay['hole_centers_from_socket'].items():
        xv = holes_xy[hid][0] - datum
        hole_err[hid] = xv - tgt
        if abs(hole_err[hid]) > tol:
            report['violations'].append(
                '%s exterior %.1f vs low A target %.1f (err %+.1f, tol %.0f)'
                % (hid, xv, tgt, hole_err[hid], tol))
    body_len_vis = body_len - datum
    len_err = body_len_vis - lay['body_length']
    if abs(len_err) > lay['body_length_tolerance']:
        report['violations'].append(
            'visible body length %.1f vs low A %.1f (err %+.1f, tol %.0f)'
            % (body_len_vis, lay['body_length'], len_err,
               lay['body_length_tolerance']))
    # reachability invariant: straight-line distance from socket to target
    # must not exceed available acoustic arc
    for hid, tgt in lay['hole_centers_from_socket'].items():
        avail = positions[hid] - body0
        if (tgt + datum) > avail + tol:
            report['violations'].append('%s low A target unreachable' % hid)

    report.update({
        'R': R, 's_bend': s_bend, 'phi6_deg': math.degrees(phi6),
        'holes_xy_body': holes_xy, 'apex_x': apex_x, 'foot_x': foot_x,
        'body_len_body': body_len, 'body_len_visible': body_len_vis,
        'len_err_vs_lowA': len_err, 'hole_err_vs_lowA': hole_err,
        'pass_separation': 2.0 * R,
        'segments': [
            {'type': 'straight', 'length': s_bend - body0},
            {'type': 'bend', 'angle_deg': 180.0, 'radius': R},
            {'type': 'straight', 'length': L_end - s_bend - math.pi * R}]})
    return report

def fold_report_str(rep):
    lines = ['fold: foot U-bend  R=%.1f mm  bend start s=%.1f  h6 %.1f deg into bend'
             % (rep['R'], rep['s_bend'], rep['phi6_deg'])]
    lines.append('  pass separation %.1f mm; apex x=%.1f; foot exit x=%.1f (return pass)'
                 % (rep['pass_separation'], rep['apex_x'], rep['foot_x']))
    lines.append('  visible body length %.1f mm (low A target %.1f, err %+.1f)'
                 % (rep['body_len_visible'],
                    rep['body_len_visible'] - rep['len_err_vs_lowA'],
                    rep['len_err_vs_lowA']))
    he = rep['hole_err_vs_lowA']
    lines.append('  hole exterior err vs low A (mm): '
                 + '  '.join('%s %+.0f' % (h, he[h])
                             for h in ['h1', 'h2', 'h3', 'h4', 'h5', 'h6']))
    lines.append('  feasible per spec tolerances: %s' % rep['feasible'])
    for v in rep['violations']:
        lines.append('    VIOLATION: ' + v)
    return '\n'.join(lines)

fold0 = solve_foot_fold(spec_straight, sc1['positions'], sc1['L_end'])
print(fold_report_str(fold0))

H6_MAX_BEND_ANGLE_DEG = 50.0   # 33*(1-cos50) ~ 12 mm lateral: offset-hole ergonomics

def solve_foot_fold(spec, positions, L_end):
    d_bore = 2.0 * bore_radius(spec)
    R = BEND_RADIUS_MIN_FACTOR * d_bore   # phi6 grows with R: legal minimum is optimal
    phi_max = math.radians(H6_MAX_BEND_ANGLE_DEG)
    s5, s6 = positions['h5'], positions['h6']
    body0 = spec['bore']['body_start_arc']
    report = {'feasible': True, 'violations': []}

    s_bend = L_end - FOOT_TAIL_MM - math.pi * R
    phi6 = max(0.0, (s6 - s_bend) / R)
    if phi6 > phi_max:
        report['violations'].append(
            'h6 %.1f deg into bend exceeds %.0f deg limit at minimum bend radius'
            % (math.degrees(phi6), H6_MAX_BEND_ANGLE_DEG))
        report['feasible'] = False
    if s_bend < s5 + H5_BEND_CLEARANCE_MM:
        report['violations'].append('bend start intrudes on h5')
        report['feasible'] = False

    holes_xy = {}
    for hid in positions:
        s = positions[hid]
        if s <= s_bend:
            holes_xy[hid] = (s - body0, 0.0)
        else:
            ph = (s - s_bend) / R
            holes_xy[hid] = ((s_bend - body0) + R * math.sin(ph),
                             R * (1.0 - math.cos(ph)))
    apex_x = (s_bend - body0) + R
    body_len = apex_x + BODY_OD / 2.0 + FOOT_CAP_CLEAR
    foot_x = (s_bend - body0) - FOOT_TAIL_MM

    lay = spec['physical_layout']
    datum = (spec['headjoint']['sleeve_end_from_window']
             - spec['headjoint']['window_to_body_socket'])
    tol = lay['hole_position_tolerance']
    hole_err = {}
    for hid, tgt in lay['hole_centers_from_socket'].items():
        xv = holes_xy[hid][0] - datum
        hole_err[hid] = xv - tgt
        if abs(hole_err[hid]) > tol:
            report['violations'].append(
                '%s exterior %.1f vs low A target %.1f (err %+.1f, tol %.0f)'
                % (hid, xv, tgt, hole_err[hid], tol))
    body_len_vis = body_len - datum
    len_err = body_len_vis - lay['body_length']
    if abs(len_err) > lay['body_length_tolerance']:
        report['violations'].append(
            'visible body length %.1f vs low A %.1f (err %+.1f, tol %.0f)'
            % (body_len_vis, lay['body_length'], len_err,
               lay['body_length_tolerance']))
    for hid, tgt in lay['hole_centers_from_socket'].items():
        avail = positions[hid] - body0
        if (tgt + datum) > avail + tol:
            report['violations'].append('%s low A target unreachable' % hid)

    report.update({
        'R': R, 's_bend': s_bend, 'phi6_deg': math.degrees(phi6),
        'holes_xy_body': holes_xy, 'apex_x': apex_x, 'foot_x': foot_x,
        'body_len_body': body_len, 'body_len_visible': body_len_vis,
        'len_err_vs_lowA': len_err, 'hole_err_vs_lowA': hole_err,
        'pass_separation': 2.0 * R,
        'segments': [
            {'type': 'straight', 'length': s_bend - body0},
            {'type': 'bend', 'angle_deg': 180.0, 'radius': R},
            {'type': 'straight', 'length': L_end - s_bend - math.pi * R}]})
    return report

def make_folded_spec(spec, positions, L_end, diameters, fold):
    out = writeback_solution(spec, positions, L_end, diameters)
    out['centerline'] = {'description':
                         'foot U-fold, h6 %.1f deg into bend' % fold['phi6_deg'],
                         'segments': fold['segments']}
    return load_spec(out)

fold1 = solve_foot_fold(spec_straight, sc1['positions'], sc1['L_end'])
print(fold_report_str(fold1))
spec_m2 = make_folded_spec(spec_straight, sc1['positions'], sc1['L_end'], S1, fold1)
print()
print('folded spec valid; re-solving with bend corrections on (kappa=1)...')
m2_pos = dict(sc1['positions'])
m2_L = sc1['L_end']
m2_pos, m2_L = optimize_sweep(spec_m2, m2_pos, m2_L, diameters=S1,
                              apply_bends=True, tol=0.4)
resid = r1_residuals(spec_m2, m2_pos, m2_L, diameters=S1, apply_bends=True)
print('after bend-aware sweep: L_end = %.2f (was %.2f)' % (m2_L, sc1['L_end']))
print('R1 residuals:', {k: round(v, 1) for k, v in resid.items()})

# iterate fold <-> bend-aware acoustics to a fixed point
for it in range(3):
    fold2 = solve_foot_fold(spec_m2, m2_pos, m2_L)
    spec_m2 = make_folded_spec(spec_m2, m2_pos, m2_L, S1, fold2)
    m2_pos, m2_L = optimize_sweep(spec_m2, m2_pos, m2_L, diameters=S1,
                                  apply_bends=True, tol=0.3)
    resid = r1_residuals(spec_m2, m2_pos, m2_L, diameters=S1, apply_bends=True)
    worst = max(abs(v) for v in resid.values())
    print('iter %d: L_end=%.2f s_bend=%.2f phi6=%.1f worst_R1=%.2fc'
          % (it, m2_L, fold2['s_bend'], fold2['phi6_deg'], worst))
    if worst < 0.5 and abs(bore_end(spec_m2) - m2_L) < 0.3:
        break

fold_final = solve_foot_fold(spec_m2, m2_pos, m2_L)
spec_m2 = make_folded_spec(spec_m2, m2_pos, m2_L, S1, fold_final)

# M1 final: straight body, S1 diameters
spec_m1 = writeback_solution(spec_straight, sc1['positions'], sc1['L_end'], S1)
res_m1 = analyze(spec_m1, apply_bends=False)
res_m2 = analyze(spec_m2, apply_bends=True)
order = [f['note'] for f in spec_m1['fingerings']]
print('M1 straight (S1 diameters):')
print(cents_table(res_m1, order))
print()
print('M2 folded (bend-corrected):')
print(cents_table(res_m2, order))
print()
print('M2 - M1 per fingering (cents):')
print('  ' + '  '.join('%s %+.1f' % (n, res_m2[n]['cents'] - res_m1[n]['cents'])
                       for n in order))
print()
print('M2 solved positions:', {k: round(v, 1) for k, v in m2_pos.items()})
print('M1 solved positions:', {k: round(v, 1) for k, v in sc1['positions'].items()})
print()
print(fold_report_str(fold_final))

solutions = {
    'm1_straight': {
        'diameters': S1,
        'positions': {k: round(v, 3) for k, v in sc1['positions'].items()},
        'L_end': round(sc1['L_end'], 3),
        'cents': {n: round(res_m1[n]['cents'], 2) for n in res_m1},
    },
    'm2_folded': {
        'diameters': S1,
        'positions': {k: round(v, 3) for k, v in m2_pos.items()},
        'L_end': round(m2_L, 3),
        'cents': {n: round(res_m2[n]['cents'], 2) for n in res_m2},
        'fold': {k: fold_final[k] for k in
                 ['R', 's_bend', 'phi6_deg', 'apex_x', 'foot_x',
                  'body_len_visible', 'len_err_vs_lowA', 'hole_err_vs_lowA',
                  'pass_separation', 'feasible', 'violations']},
    },
    'scenarios': {
        'S0_seed': 'worst R2 C#6 -97c, span 209mm',
        'S1_big': 'worst R2 C#6 -26c, span 205mm  <- selected',
        'S2_compress': 'worst R2 F#5 -37c, span 180mm',
        'S3_mid': 'worst R2 C#6 -24c, span 194mm',
    },
    'bend_correction': 'toroidal inertance ratio (1+sqrt(1-B^2))/2, kappa=1.0, '
                       'Nederveen JASA 104:1616 (1998); Felix/Dalmont/Nederveen '
                       'JASA 131:4164 (2012)',
}
p = save_json('foldflute_solutions', solutions)
print('saved:', p)

from build123d import *
# ---------------- derived CAD parameters (mm) from solver state ----------------
BODY0    = spec_m2['bore']['body_start_arc']          # 120.0 arc datum
BORE_R   = bore_radius(spec_m2)                       # 11.0
BODY_ARC = m2_L - BODY0                               # total body air path
SBEND    = fold_final['s_bend'] - BODY0               # bend start, body x
BR       = fold_final['R']                            # 33.0 bend radius
TAIL     = m2_L - fold_final['s_bend'] - math.pi * BR # 6.0 return tail
OD_R     = BODY_OD / 2.0                              # 14.3 body outer radius
SPIG_R   = spec_m2['headjoint']['printed_body_mating_outer_diameter'] / 2.0  # 12.4
SPIG_LEN = spec_m2['headjoint']['sleeve_extension_beyond_tube_end']   # 56.0
JOINT_X  = 150.0        # A/B barrel joint (no holes above; B fits print bed)
TENON_R, TENON_LEN = 13.0, 20.0        # A male tenon OD 26.0
SOCK_R, SOCK_DEPTH = 13.1, 20.2        # B socket ID 26.2 (0.2 diametral fit)
RING_R, RING_END = 15.4, 180.0         # B reinforcing ring OD 30.8 over socket
FOOT_X   = SBEND - TAIL                # foot opening x on return pass
PASS_Y   = 2.0 * BR                    # 66.0 return pass centreline y

hole_geo = {}   # id -> (x, y, r_hole)
for h in spec_m2['holes']:
    s = h['position'] - BODY0
    if s <= SBEND:
        hx, hy = s, 0.0
    else:
        ph = (s - SBEND) / BR
        hx, hy = SBEND + BR * math.sin(ph), BR * (1.0 - math.cos(ph))
    hole_geo[h['id']] = (hx, hy, h['diameter'] / 2.0)
print('BODY_ARC %.2f  SBEND %.2f  TAIL %.2f  FOOT_X %.2f' %
      (BODY_ARC, SBEND, TAIL, FOOT_X))
print('holes:', {k: (round(v[0], 1), round(v[1], 1), v[2]) for k, v in hole_geo.items()})

# ---------------- canonical air column ----------------
with BuildLine() as bore_path:
    l1 = Line((0, 0, 0), (SBEND, 0, 0))
    a1 = JernArc(start=(SBEND, 0), tangent=(1, 0), radius=BR, arc_size=180)
    l2 = Line((SBEND, PASS_Y, 0), (FOOT_X, PASS_Y, 0))
with BuildSketch(Plane.YZ) as bore_prof:
    Circle(BORE_R)
air = sweep(bore_prof.sketch, path=bore_path.line)
show(air, 'air_solid')

# build123d's star-import clobbered MM (unit constant). Re-bind the adapter's
# mm->m factor to a private name so CAD and acoustics can't collide again.
MM_TO_M = 0.001

def to_openwind_geometry(spec, positions=None, diameters=None, chimneys=None,
                         L_end=None, apply_bends=True, kappa=BEND_KAPPA_DEFAULT,
                         extra_holes=None):
    cal = spec['calibration']
    lead = spec['headjoint']['effective_length_correction']
    end = L_end if L_end is not None else bore_end(spec)
    smap = arc_effective_map(spec, kappa) if apply_bends else (lambda s: s)
    end_eff = smap(end) + cal['global_length_correction'] + lead
    r_m = bore_radius(spec) * MM_TO_M
    main_bore = [[0.0, r_m], [end_eff * MM_TO_M, r_m]]
    rows = [['label', 'position', 'radius', 'chimney']]
    hole_states = {}
    for h in effective_holes(spec, positions, diameters, chimneys):
        hid = h['id']
        rad = h['diameter'] / 2.0 + cal['hole_radius_corrections'].get(hid, 0.0)
        chy = h['chimney'] + cal['chimney_corrections'].get(hid, 0.0)
        pos = smap(h['position']) + lead
        rows.append([hid, pos * MM_TO_M, rad * MM_TO_M, chy * MM_TO_M])
        hole_states[hid] = None
    for h in (extra_holes or []):
        pos = smap(h['position']) + lead
        rows.append([h['id'], pos * MM_TO_M, h['diameter'] / 2.0 * MM_TO_M,
                     h['chimney'] * MM_TO_M])
        hole_states[h['id']] = h['state']
    notes = [ow_label(f['note']) for f in spec['fingerings']]
    chart = [['label'] + notes]
    state_char = {'closed': 'x', 'open': 'o'}
    for hid in hole_states:
        if hole_states[hid] is None:
            row = [hid]
            for f in spec['fingerings']:
                row.append(state_char[f['holes'][hid]])
        else:
            row = [hid] + [state_char[hole_states[hid]]] * len(notes)
        chart.append(row)
    return main_bore, rows, chart, cal['temperature_c']

# regression check: analyze must reproduce the M2 D4 solution
chk = analyze(spec_m2, notes=['D4'], apply_bends=True)
print('adapter OK, D4 cents = %+.2f (expect ~0)' % chk['D4']['cents'])

with BuildLine() as body_path:
    Line((SPIG_LEN, 0, 0), (SBEND, 0, 0))
    JernArc(start=(SBEND, 0), tangent=(1, 0), radius=BR, arc_size=180)
    Line((SBEND, PASS_Y, 0), (FOOT_X, PASS_Y, 0))
with BuildSketch(Plane.YZ.offset(SPIG_LEN)) as body_prof:
    Circle(OD_R)
main_tube = sweep(body_prof.sketch, path=body_path.line)

spigot = Cylinder(SPIG_R, SPIG_LEN, rotation=(0, 90, 0)).move(
    Location((SPIG_LEN / 2.0, 0, 0)))
ring = Cylinder(RING_R, RING_END - JOINT_X, rotation=(0, 90, 0)).move(
    Location(((JOINT_X + RING_END) / 2.0, 0, 0)))
body_raw = spigot.fuse(main_tube).fuse(ring)
show(body_raw, 'body_raw')

body = body_raw - air
show(body, 'body_solid')

cutters = None
for hid in ['h1', 'h2', 'h3', 'h4', 'h5', 'h6']:
    hx, hy, hr = hole_geo[hid]
    c = Cylinder(hr, 14.0, rotation=(0, 0, 0)).move(Location((hx, hy, 6.0 + 7.0)))
    cutters = c if cutters is None else cutters.fuse(c)
body_holed = body - cutters
show(body_holed, 'body_solid')

JOINT_X = 172.0            # between h1 (edge 169.7) and h2 (edge 194.7)
TENON_LEN = 20.0
SOCK_DEPTH = 20.4          # 0.4 axial slack
RING_X0, RING_X1 = 172.0, 194.0

ring = Cylinder(RING_R, RING_X1 - RING_X0, rotation=(0, 90, 0)).move(
    Location(((RING_X0 + RING_X1) / 2.0, 0, 0)))
body_raw = spigot.fuse(main_tube).fuse(ring)
body = body_raw - air
cutters = None
for hid in ['h1', 'h2', 'h3', 'h4', 'h5', 'h6']:
    hx, hy, hr = hole_geo[hid]
    c = Cylinder(hr, 14.0).move(Location((hx, hy, 13.0)))
    cutters = c if cutters is None else cutters.fuse(c)
body = body - cutters
show(body, 'body_solid')

CHIMNEY_MEASURED = OD_R - BORE_R    # 3.30 at hole axis, from as-built CAD
chims = {h['id']: CHIMNEY_MEASURED for h in spec_m2['holes']}

pos_before = dict(m2_pos)
m2_pos, m2_L = optimize_sweep(spec_m2, m2_pos, m2_L, diameters=S1,
                              chimneys=chims, apply_bends=True, tol=0.3)
resid = r1_residuals(spec_m2, m2_pos, m2_L, diameters=S1, chimneys=chims,
                     apply_bends=True)
worst = max(abs(v) for v in resid.values())
drift = {k: round(m2_pos[k] - pos_before[k], 2) for k in m2_pos}
print('chimney 3.0 -> 3.3 re-solve: worst R1 %.2f c' % worst)
print('position drift (mm):', drift, ' L_end drift: %+.2f' % (m2_L - bore_end(spec_m2)))

fold_final = solve_foot_fold(spec_m2, m2_pos, m2_L)
spec_m2 = make_folded_spec(spec_m2, m2_pos, m2_L, S1, fold_final)
for h in spec_m2['holes']:
    h['chimney_height_measured'] = round(CHIMNEY_MEASURED, 3)

res_m2 = analyze(spec_m2, apply_bends=True)
worst_all = max(abs(res_m2[n]['cents'])
                for n in ['D4', 'E4', 'F#4', 'G4', 'A4', 'B4', 'C#5'])
print('final M2: L_end %.2f  s_bend %.2f  phi6 %.1f  worst R1 %.2f c'
      % (m2_L, fold_final['s_bend'], fold_final['phi6_deg'], worst_all))

# ---- final parametric rebuild at converged coordinates ----
BODY_ARC = m2_L - BODY0
SBEND = fold_final['s_bend'] - BODY0
TAIL = m2_L - fold_final['s_bend'] - math.pi * BR
FOOT_X = SBEND - TAIL
hole_geo = {}
for h in spec_m2['holes']:
    s = h['position'] - BODY0
    if s <= SBEND:
        hx, hy = s, 0.0
    else:
        ph = (s - SBEND) / BR
        hx, hy = SBEND + BR * math.sin(ph), BR * (1.0 - math.cos(ph))
    hole_geo[h['id']] = (hx, hy, h['diameter'] / 2.0)

with BuildLine() as bore_path:
    Line((0, 0, 0), (SBEND, 0, 0))
    JernArc(start=(SBEND, 0), tangent=(1, 0), radius=BR, arc_size=180)
    Line((SBEND, PASS_Y, 0), (FOOT_X, PASS_Y, 0))
with BuildSketch(Plane.YZ) as bore_prof:
    Circle(BORE_R)
air = sweep(bore_prof.sketch, path=bore_path.line)
show(air, 'air_solid')

with BuildLine() as body_path:
    Line((SPIG_LEN, 0, 0), (SBEND, 0, 0))
    JernArc(start=(SBEND, 0), tangent=(1, 0), radius=BR, arc_size=180)
    Line((SBEND, PASS_Y, 0), (FOOT_X, PASS_Y, 0))
with BuildSketch(Plane.YZ.offset(SPIG_LEN)) as body_prof:
    Circle(OD_R)
main_tube = sweep(body_prof.sketch, path=body_path.line)
ring = Cylinder(RING_R, RING_X1 - RING_X0, rotation=(0, 90, 0)).move(
    Location(((RING_X0 + RING_X1) / 2.0, 0, 0)))
body_raw = spigot.fuse(main_tube).fuse(ring)
body = body_raw - air
cutters = None
for hid in ['h1', 'h2', 'h3', 'h4', 'h5', 'h6']:
    hx, hy, hr = hole_geo[hid]
    c = Cylinder(hr, 14.0).move(Location((hx, hy, 13.0)))
    cutters = c if cutters is None else cutters.fuse(c)
body = body - cutters
show(body, 'body_solid')
print('holes final:', {k: (round(v[0], 2), round(v[1], 2)) for k, v in hole_geo.items()})

# seam alignment bosses on the canonical body (single canonical solid first)
BOSS_R, BOSS_H = 4.5, 5.2
boss_xy = [(215.0, 16.0), (300.0, -16.0), (394.2, 33.0)]
for bx, by in boss_xy:
    body = body.fuse(Cylinder(BOSS_R, BOSS_H).move(Location((bx, by, 0))))
show(body, 'body_solid')

# ---- section A: monolithic upper barrel with male tenon ----
BIG = 1000.0
cutA = Box(BIG, BIG, BIG, align=(Align.MAX, Align.CENTER, Align.CENTER)).move(
    Location((JOINT_X, 0, 0)))
tenon = (Cylinder(TENON_R, TENON_LEN, rotation=(0, 90, 0))
         - Cylinder(BORE_R, TENON_LEN, rotation=(0, 90, 0))).move(
    Location((JOINT_X + TENON_LEN / 2.0, 0, 0)))
sectionA = (body & cutA).fuse(tenon)
show(sectionA, 'section_a')

cutB = Box(BIG, BIG, BIG, align=(Align.MIN, Align.CENTER, Align.CENTER)).move(
    Location((JOINT_X, 0, 0)))
socket_recess = Cylinder(SOCK_R, SOCK_DEPTH, rotation=(0, 90, 0)).move(
    Location((JOINT_X + SOCK_DEPTH / 2.0, 0, 0)))
sectionB = (body & cutB) - socket_recess
show(sectionB, 'section_b')

# assembled-fit check data, then shells
c = sectionA.bounding_box()
print('A bbox x: %.1f..%.1f' % (c.min.X, c.max.X))
c = sectionB.bounding_box()
print('B bbox x: %.1f..%.1f  y: %.1f..%.1f  z: %.1f..%.1f'
      % (c.min.X, c.max.X, c.min.Y, c.max.Y, c.min.Z, c.max.Z))

halfF = Box(BIG, BIG, BIG, align=(Align.CENTER, Align.CENTER, Align.MIN))
halfK = Box(BIG, BIG, BIG, align=(Align.CENTER, Align.CENTER, Align.MAX))
PIN_HOLE_R, PIN_HOLE_DEPTH = 1.25, 2.2
shell_front = sectionB & halfF
shell_back = sectionB & halfK
for bx, by in boss_xy:
    shell_front = shell_front - Cylinder(PIN_HOLE_R, PIN_HOLE_DEPTH).move(
        Location((bx, by, PIN_HOLE_DEPTH / 2.0)))
    shell_back = shell_back - Cylinder(PIN_HOLE_R, PIN_HOLE_DEPTH).move(
        Location((bx, by, -PIN_HOLE_DEPTH / 2.0)))
show(shell_front, 'left_shell')    # +z half: carries all section-B toneholes
show(shell_back, 'right_shell')
pin = Cylinder(1.15, 4.0)
show(pin, 'alignment_pin')

# bed-oriented copies (z=0 = print bed)
a_print = sectionA.rotate(Axis.Y, 90).move(Location((0, 0, 192.0)))
bb = a_print.bounding_box()
print('A print bbox z: %.1f..%.1f xy: %.0fx%.0f'
      % (bb.min.Z, bb.max.Z, bb.size.X, bb.size.Y))
show(a_print, 'section_a_print')

back_print = shell_back.rotate(Axis.X, 180)
bb = back_print.bounding_box()
print('back shell print bbox z: %.1f..%.1f xy: %.0fx%.0f'
      % (bb.min.Z, bb.max.Z, bb.size.X, bb.size.Y))
show(back_print, 'right_shell_print')

bb = shell_front.bounding_box()
print('front shell print bbox z: %.1f..%.1f xy: %.0fx%.0f'
      % (bb.min.Z, bb.max.Z, bb.size.X, bb.size.Y))

pin_print = pin.move(Location((0, 0, 2.0)))
show(pin_print, 'pin_print')

# --- M1 straight: re-solve with measured chimney 3.3 ---
m1_pos = dict(sc1['positions'])
m1_Lf = sc1['L_end']
m1_pos, m1_Lf = optimize_sweep(spec_straight, m1_pos, m1_Lf, diameters=S1,
                               chimneys=chims, apply_bends=False, tol=0.3)
resid = r1_residuals(spec_straight, m1_pos, m1_Lf, diameters=S1,
                     chimneys=chims, apply_bends=False)
print('M1 chimney-corrected: L_end %.2f worst R1 %.2f c'
      % (m1_Lf, max(abs(v) for v in resid.values())))
spec_m1 = writeback_solution(spec_straight, m1_pos, m1_Lf, S1)
for h in spec_m1['holes']:
    h['chimney_height_measured'] = round(CHIMNEY_MEASURED, 3)
res_m1 = analyze(spec_m1, apply_bends=False)
print('M1 holes body-x:', {k: round(v - BODY0, 1) for k, v in m1_pos.items()})

# --- M1 straight body CAD: two vertical-print barrel sections ---
M1_ARC = m1_Lf - BODY0                 # 451.60 body length (bore open at foot)
M1_JOINT = 219.0                       # between h2 (edge 204.5) and h3 (edge 233.9)
M1_RING0, M1_RING1 = 219.0, 241.0

tube_m1 = Cylinder(OD_R, M1_ARC - SPIG_LEN, rotation=(0, 90, 0)).move(
    Location(((SPIG_LEN + M1_ARC) / 2.0, 0, 0)))
ring_m1 = Cylinder(RING_R, M1_RING1 - M1_RING0, rotation=(0, 90, 0)).move(
    Location(((M1_RING0 + M1_RING1) / 2.0, 0, 0)))
air_m1 = Cylinder(BORE_R, M1_ARC, rotation=(0, 90, 0)).move(
    Location((M1_ARC / 2.0, 0, 0)))
body_m1 = spigot.fuse(tube_m1).fuse(ring_m1) - air_m1
cutters = None
for hid in ['h1', 'h2', 'h3', 'h4', 'h5', 'h6']:
    c = Cylinder(S1[hid] / 2.0, 14.0).move(
        Location((m1_pos[hid] - BODY0, 0, 13.0)))
    cutters = c if cutters is None else cutters.fuse(c)
body_m1 = body_m1 - cutters
show(body_m1, 'm1_body')

cutA1 = Box(BIG, BIG, BIG, align=(Align.MAX, Align.CENTER, Align.CENTER)).move(
    Location((M1_JOINT, 0, 0)))
cutB1 = Box(BIG, BIG, BIG, align=(Align.MIN, Align.CENTER, Align.CENTER)).move(
    Location((M1_JOINT, 0, 0)))
tenon1 = (Cylinder(TENON_R, TENON_LEN, rotation=(0, 90, 0))
          - Cylinder(BORE_R, TENON_LEN, rotation=(0, 90, 0))).move(
    Location((M1_JOINT + TENON_LEN / 2.0, 0, 0)))
sock1 = Cylinder(SOCK_R, SOCK_DEPTH, rotation=(0, 90, 0)).move(
    Location((M1_JOINT + SOCK_DEPTH / 2.0, 0, 0)))
m1_a = (body_m1 & cutA1).fuse(tenon1)
m1_b = (body_m1 & cutB1) - sock1
show(m1_a, 'm1_section_a')
show(m1_b, 'm1_section_b')
print('m1_a bbox x 0..%.1f   m1_b x %.1f..%.1f'
      % (m1_a.bounding_box().max.X, m1_b.bounding_box().min.X,
         m1_b.bounding_box().max.X))

m1a_print = m1_a.rotate(Axis.Y, 90).move(Location((0, 0, 239.0)))
show(m1a_print, 'm1_a_print')
m1b_print = m1_b.rotate(Axis.Y, -90).move(Location((0, 0, -219.0)))
bb = m1b_print.bounding_box()
print('m1_b print z %.1f..%.1f' % (bb.min.Z, bb.max.Z))
show(m1b_print, 'm1_b_print')

order = [f['note'] for f in spec_m1['fingerings']]
final = {
    'm1_straight': {
        'spec': spec_m1,
        'cents': {n: round(res_m1[n]['cents'], 2) for n in order},
        'sections': {'a': '0-219 body-x +20 tenon, vertical print, holes h1-h2',
                     'b': '219-451.6 body-x, vertical print, holes h3-h6'},
    },
    'm2_folded': {
        'spec': spec_m2,
        'cents': {n: round(res_m2[n]['cents'], 2) for n in order},
        'fold': {k: fold_final[k] for k in
                 ['R', 's_bend', 'phi6_deg', 'apex_x', 'foot_x',
                  'body_len_visible', 'len_err_vs_lowA', 'hole_err_vs_lowA',
                  'pass_separation', 'feasible', 'violations']},
        'sections': {'a': '0-172 body-x +20 tenon, vertical print, hole h1',
                     'b': '172-end incl U-lobe, split shells, holes h2-h6'},
    },
    'chimney_measured_mm': round(CHIMNEY_MEASURED, 3),
}
print(save_json('foldflute_final', final))
m2_minus_m1 = {n: round(res_m2[n]['cents'] - res_m1[n]['cents'], 2) for n in order}
print('M2-M1 spread:', min(m2_minus_m1.values()), 'to', max(m2_minus_m1.values()))