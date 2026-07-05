# foldflute prelude -- deterministic acoustic + fold-path helpers (SPEC.md 6-7).
#
# This file is designed to run INSIDE the build123d-mcp execute() sandbox
# (no file/network access; imports restricted to the server allowlist), so:
#   - the InstrumentSpec arrives as a plain dict (load it host-side from the
#     JSONC file and inject it into the session);
#   - no dataclasses, no open(), no getattr(), no explicit dunder access;
#   - build123d's star-import defines MM=1, so the mm->m factor here is named
#     MM_TO_M to avoid collision.
#
# Units: mm and Hz at every public boundary. SI (m) only inside
# to_openwind_geometry, at the openwind boundary (SPEC.md section 5).

import math
import copy
import numpy as np
from openwind import ImpedanceComputation

# ---------------------------------------------------------------- notes
NOTE_SEMITONES = {'C': -9, 'C#': -8, 'Db': -8, 'D': -7, 'D#': -6, 'Eb': -6,
                  'E': -5, 'F': -4, 'F#': -3, 'Gb': -3, 'G': -2, 'G#': -1,
                  'Ab': -1, 'A': 0, 'A#': 1, 'Bb': 1, 'B': 2}

def note_hz(name, pitch_standard_hz=440.0):
    """Equal-temperament frequency of a note name like 'F#4'."""
    letter = name[:-1]
    octave = int(name[-1])
    semis = NOTE_SEMITONES[letter] + 12 * (octave - 4)
    return pitch_standard_hz * (2.0 ** (semis / 12.0))

def cents_of(f, f_ref):
    return 1200.0 * math.log(f / f_ref) / math.log(2.0)

def ow_label(note):
    """Openwind chart labels cannot contain '#' (comment character)."""
    return note.replace('#', 's')

# ---------------------------------------------------------------- spec
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
    """Validate an InstrumentSpec dict against the SPEC.md section 5 invariants.

    position_bounds may extend past the (retunable) bore end; only optimized
    hole positions are required to sit inside the bore.
    """
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
        if positions and hid in positions:
            pos = positions[hid]
        if diameters and hid in diameters:
            dia = diameters[hid]
        if chimneys and hid in chimneys:
            chy = chimneys[hid]
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

# ------------------------------------------------ bend corrections (SPEC 6.2)
# A toroidal bend of centerline radius R and bore radius a behaves, in the
# long-wavelength limit, as a straight duct of the same volume but reduced
# inertance: with B = a/R and a 1/r potential-flow velocity profile across the
# section,
#     M_bend / M_straight = (1 + sqrt(1 - B^2)) / 2  ~=  1 - B^2/4.
# Reduced inertance at unchanged compliance makes the bend acoustically
# SHORTER than its arc length (pitch rises). The true correction depends on
# the bend's position in the standing wave (inertance- vs compliance-
# dominated); kappa = 1 applies the full velocity-antinode value and is the
# calibration-refinable knob.
# Refs: C.J. Nederveen, JASA 104(3):1616-1626 (1998); S. Felix, J.-P. Dalmont,
# C.J. Nederveen, JASA 131:4164-4172 (2012).
BEND_KAPPA_DEFAULT = 1.0

def bend_shrink_factor(bore_radius_mm, bend_radius_mm, kappa=BEND_KAPPA_DEFAULT):
    """Fraction of bend arc length removed from the effective air column."""
    B = bore_radius_mm / bend_radius_mm
    if B >= 1.0:
        raise ValueError('bend radius must exceed bore radius')
    ratio = (1.0 + math.sqrt(1.0 - B * B)) / 2.0
    return kappa * (1.0 - ratio)

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

# ------------------------------------------------ openwind adapter + analyze
MM_TO_M = 0.001

def to_openwind_geometry(spec, positions=None, diameters=None, chimneys=None,
                         L_end=None, apply_bends=True, kappa=BEND_KAPPA_DEFAULT,
                         extra_holes=None):
    """InstrumentSpec (mm) -> (main_bore_m, holes_table_m, chart, temperature_c).

    Applies the headjoint effective-length correction, calibration
    corrections, and toroidal bend corrections. extra_holes: list of
    {'id','position','diameter','chimney','state'} always-open/closed extras
    (e.g. a drain-hole experiment).
    """
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
    """Frequency of the |Z| minimum nearest target_hz (whistle resonance).

    Openwind's phase-based antiresonance_frequencies() proved unreliable on
    coarse grids, so this uses a local geometric grid around the target,
    widened if the minimum lands on a boundary, with two parabolic
    refinement passes in (log f, log |Z|).
    """
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
        spec, positions, diameters, chimneys, L_end, apply_bends, kappa,
        extra_holes)
    first = chart[0][1]
    fs0 = np.array([200.0, 300.0, 400.0])
    return ImpedanceComputation(fs0, main_bore, holes, chart, note=first,
                                temperature=temp, losses=True,
                                radiation_category='unflanged')

def analyze(spec, notes=None, positions=None, diameters=None, chimneys=None,
            L_end=None, apply_bends=True, kappa=BEND_KAPPA_DEFAULT,
            extra_holes=None, comp=None):
    """Passive resonance prediction per fingering (impedance minima, SPEC 6).

    Returns {note: {'target_hz', 'f_res', 'cents'}}.
    """
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

def cents_table(res, order=None):
    keys = order if order else list(res.keys())
    lines = ['note    target(Hz)  res(Hz)   cents']
    for n in keys:
        r = res[n]
        lines.append('%-6s  %9.2f  %8.2f  %+6.1f'
                     % (n, r['target_hz'], r['f_res'], r['cents']))
    return '\n'.join(lines)

# ---------------------------------------------------------------- optimizer
# Gauss-Seidel over (bore end -> D4) then (h6..h1 -> E4..C#5), damped secant
# per scalar parameter. Register 2 and D6 are reported as predictions, not
# targets: register balance is set by the diameter choices, which are inputs.
R1_TARGET_NOTES = [('h6', 'E4'), ('h5', 'F#4'), ('h4', 'G4'),
                   ('h3', 'A4'), ('h2', 'B4'), ('h1', 'C#5')]
L_END_BOUNDS = (545.0, 608.0)
HOLE_END_MARGIN = 8.0   # mm: keep the last hole clear of the open foot

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
    """One Gauss-Seidel sweep; mutates positions, returns (positions, L_end).

    Call repeatedly until r1_residuals() is inside tolerance (2 sweeps
    typically reach +/-0.3 cents from seed positions).
    """
    fnL = lambda L: eval_note_cents(spec, 'D4', positions, L, diameters,
                                    chimneys, apply_bends, kappa)
    L_end, _ = _secant_solve(fnL, L_end, L_END_BOUNDS, slope0=-3.0, tol=tol)
    hole_bounds = {h['id']: (h['position_bounds'][0],
                             min(h['position_bounds'][1],
                                 L_end - HOLE_END_MARGIN))
                   for h in spec['holes']}
    for hid, note in R1_TARGET_NOTES:
        fn = (lambda p, hid=hid, note=note:
              eval_note_cents(spec, note, {**positions, hid: p}, L_end,
                              diameters, chimneys, apply_bends, kappa))
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

# ---------------------------------------------------------------- fold path
# Topology search result for the 22 mm cylindrical low D bore (see
# reports/REPORT.md): a 180-degree bend needs pi*R >= 104 mm of arc at the
# minimum legal radius (1.5 x bore diameter), which fits neither between
# adjacent toneholes (< 57 mm gaps), nor above h1 (~80 mm of straight
# available vs >= 207 mm needed for an order-preserving fold-back), nor
# wholly below h6 (~90 mm of tail). The one legal fold is a foot U-bend with
# h6 a bounded angle into the bend. The fold plane is XY (passes side by
# side laterally); every hole exit stays +z as the spec exit normals require.
BEND_RADIUS_MIN_FACTOR = 1.5
FOOT_TAIL_MM = 6.0          # bore continuing past the bend to the foot opening
H6_MAX_BEND_ANGLE_DEG = 50.0   # R*(1-cos 50) ~ 12 mm lateral: offset-hole ergonomics
H5_BEND_CLEARANCE_MM = 12.0
BODY_OD = 28.6              # body outer diameter; matches the sleeve OD
FOOT_CAP_CLEAR = 2.0

def solve_foot_fold(spec, positions, L_end):
    """Solve the foot U-fold. Hard constraints (bend radius, h6 bend angle,
    h5 clearance) set 'feasible'; low A layout misses are reported in
    'layout_violations' and handled by the spec's relax_visual_length policy.
    """
    d_bore = 2.0 * bore_radius(spec)
    R = BEND_RADIUS_MIN_FACTOR * d_bore   # phi6 grows with R: legal minimum is optimal
    phi_max = math.radians(H6_MAX_BEND_ANGLE_DEG)
    s5, s6 = positions['h5'], positions['h6']
    body0 = spec['bore']['body_start_arc']
    report = {'feasible': True, 'violations': [], 'layout_violations': []}

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
    foot_x = (s_bend - body0) - FOOT_TAIL_MM   # on the return pass, y = 2R

    # Low A layout comparison. Datum assumption: hole_centers_from_socket are
    # measured from where the body emerges from the sleeve (sleeve end), i.e.
    # visible x = body x - (sleeve_end_from_window - window_to_body_socket).
    lay = spec['physical_layout']
    datum = (spec['headjoint']['sleeve_end_from_window']
             - spec['headjoint']['window_to_body_socket'])
    tol = lay['hole_position_tolerance']
    hole_err = {}
    for hid, tgt in lay['hole_centers_from_socket'].items():
        xv = holes_xy[hid][0] - datum
        hole_err[hid] = xv - tgt
        if abs(hole_err[hid]) > tol:
            report['layout_violations'].append(
                '%s exterior %.1f vs low A target %.1f (err %+.1f, tol %.0f)'
                % (hid, xv, tgt, hole_err[hid], tol))
    body_len_vis = body_len - datum
    len_err = body_len_vis - lay['body_length']
    if abs(len_err) > lay['body_length_tolerance']:
        report['layout_violations'].append(
            'visible body length %.1f vs low A %.1f (err %+.1f, tol %.0f)'
            % (body_len_vis, lay['body_length'], len_err,
               lay['body_length_tolerance']))
    for hid, tgt in lay['hole_centers_from_socket'].items():
        avail = positions[hid] - body0
        if (tgt + datum) > avail + tol:
            report['layout_violations'].append(
                '%s low A target unreachable' % hid)

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

# --------------------------------------------- miter elbows (Coltman 2006)
# For internal/convoluted folds where toroidal bends cannot fit. A 90-deg
# mitered bend with its outer corner beveled at 45 deg (bevel opening
# d = 1.26*ID) has its characteristic impedance restored, so it behaves as a
# straight tube SHORTER than its centerline by a constant, mode-independent
# 0.32*ID per elbow. Refs: J.W. Coltman, "Acoustic properties of miter bends"
# (2006) eq. 3-4; Dequand et al., Acta Acustica 89:1025 (2003).
MITER_SHORTEN_FACTOR = 0.32
MITER_BEVEL_D_FACTOR = 1.26

def miter_elbow_shorten(spec):
    return MITER_SHORTEN_FACTOR * 2.0 * bore_radius(spec)

def elbow_effective_map(elbow_arcs, shorten_each):
    """Physical arc (mm from window) -> effective arc with a point shortening
    at each compensated miter elbow."""
    arcs = sorted(elbow_arcs)
    def f(s):
        return s - shorten_each * sum(1 for e in arcs if e < s)
    return f

def to_openwind_geometry_v2(spec, positions, chimneys, L_end, elbow_arcs,
                            diameters=None):
    """As to_openwind_geometry, but positions/L_end are PHYSICAL arcs on a
    convoluted (mitered) bore path; elbow corrections map them to effective.
    Used by the internal-fold feasibility study (reports/
    convolution_feasibility.json): graded deep chimneys can place holes at
    low A exterior targets, but register 2 separates by up to 1.5 semitones —
    see REPORT.md section 9 before reusing this path for a build."""
    cal = spec['calibration']
    lead = spec['headjoint']['effective_length_correction']
    smap = elbow_effective_map(elbow_arcs, miter_elbow_shorten(spec))
    end_eff = smap(L_end) + cal['global_length_correction'] + lead
    r_m = bore_radius(spec) * MM_TO_M
    main_bore = [[0.0, r_m], [end_eff * MM_TO_M, r_m]]
    rows = [['label', 'position', 'radius', 'chimney']]
    for h in effective_holes(spec, positions, diameters, chimneys):
        rows.append([h['id'], (smap(h['position']) + lead) * MM_TO_M,
                     h['diameter'] / 2.0 * MM_TO_M, h['chimney'] * MM_TO_M])
    notes = [ow_label(f['note']) for f in spec['fingerings']]
    chart = [['label'] + notes]
    state_char = {'closed': 'x', 'open': 'o'}
    for h in spec['holes']:
        row = [h['id']]
        for f in spec['fingerings']:
            row.append(state_char[f['holes'][h['id']]])
        chart.append(row)
    return main_bore, rows, chart, cal['temperature_c']

def analyze_v2(spec, positions, chimneys, L_end, elbow_arcs, diameters=None,
               notes=None):
    main_bore, holes, chart, temp = to_openwind_geometry_v2(
        spec, positions, chimneys, L_end, elbow_arcs, diameters)
    comp = ImpedanceComputation(np.array([200.0, 300.0, 400.0]), main_bore,
                                holes, chart, note=chart[0][1],
                                temperature=temp, losses=True,
                                radiation_category='unflanged')
    std = spec['target']['pitch_standard_hz']
    want = notes if notes is not None else [f['note'] for f in spec['fingerings']]
    out = {}
    for n in want:
        t = note_hz(n, std)
        comp.set_note(ow_label(n))
        out[n] = {'target_hz': t}
        f_res = _find_impedance_min(comp, t)
        out[n].update({'f_res': f_res, 'cents': cents_of(f_res, t)})
    return out

def fold_report_str(rep):
    lines = ['fold: foot U-bend  R=%.1f mm  bend start s=%.1f  h6 %.1f deg into bend'
             % (rep['R'], rep['s_bend'], rep['phi6_deg'])]
    lines.append('  pass separation %.1f mm; apex x=%.1f; foot exit x=%.1f (return pass)'
                 % (rep['pass_separation'], rep['apex_x'], rep['foot_x']))
    lines.append('  visible body length %.1f mm (err vs low A %+.1f)'
                 % (rep['body_len_visible'], rep['len_err_vs_lowA']))
    he = rep['hole_err_vs_lowA']
    lines.append('  hole exterior err vs low A (mm): '
                 + '  '.join('%s %+.0f' % (h, he[h])
                             for h in ['h1', 'h2', 'h3', 'h4', 'h5', 'h6']))
    lines.append('  hard constraints feasible: %s' % rep['feasible'])
    for v in rep['violations']:
        lines.append('    HARD: ' + v)
    for v in rep['layout_violations']:
        lines.append('    LAYOUT (relaxed per policy): ' + v)
    return '\n'.join(lines)
