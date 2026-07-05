"""foldflute geometry kernel — the single source of geometric truth.

Pure Python (math + json only; numpy optional and unused here) so it imports
unchanged in every environment: the build123d CAD venv, the NGSolve FEM venv,
and the sandboxed build123d-mcp session.

Everything geometric about a design is DERIVED here from one canonical
document (spec/*.json):

  - r_bore(x)        bore radius along the acoustic axis (mm from window)
  - face_plane()     the flat exterior face as z = z0 + slope*x
  - hole_exit(h)     where a tonehole's angled chimney meets the face
  - manifest()       named openings (window / holes / foot) for FEM face-tagging
  - to_openwind()    the 1-D transmission-line inputs for a fingering

No other module recomputes bore radius, chimney paths, or exit points. If a
number about the geometry is needed anywhere (CAD, FEM, acoustics, invariant
checks), it comes from here, so the four descriptions that used to disagree
are now one.

Axes (mm): x runs along the acoustic bore from the window datum (x=0); the
bore axis lies on z=0; the flat front face is below at z<0; the chimney tilt
angle is measured from the face normal (+z is out of the face), positive =
leaning toward the foot (down-bore).
"""
import copy
import json
import math

HOLE_IDS = ['h1', 'h2', 'h3', 'h4', 'h5', 'h6']
_STATE = {'closed': 'x', 'open': 'o'}
NOTE_SEMITONES = {'C': -9, 'C#': -8, 'Db': -8, 'D': -7, 'D#': -6, 'Eb': -6,
                  'E': -5, 'F': -4, 'F#': -3, 'Gb': -3, 'G': -2, 'G#': -1,
                  'Ab': -1, 'A': 0, 'A#': 1, 'Bb': 1, 'B': 2}


# ----------------------------------------------------------------- notes
def note_hz(name, pitch_standard_hz=440.0):
    return pitch_standard_hz * 2.0 ** (
        (NOTE_SEMITONES[name[:-1]] + 12 * (int(name[-1]) - 4)) / 12.0)


def cents(f, f_ref):
    return 1200.0 * math.log(f / f_ref, 2.0)


def ow_label(note):
    return note.replace('#', 's')          # '#' is a comment char in openwind


# ----------------------------------------------------------------- load
def load(path):
    with open(path) as fh:
        doc = json.load(fh)
    validate(doc)
    return doc


def hole(doc, hid):
    for h in doc['holes']:
        if h['id'] == hid:
            return h
    raise KeyError(hid)


# ------------------------------------------------------------- bore + face
def r_bore(doc, x):
    """Bore radius at acoustic position x (mm). Cylindrical to body_start,
    then a linear taper socket_radius -> foot_radius at the foot."""
    b = doc['bore']
    if x <= b['body_start']:
        return b['socket_radius']
    f = (x - b['body_start']) / (b['length'] - b['body_start'])
    f = min(max(f, 0.0), 1.0)
    return b['socket_radius'] + f * (b['foot_radius'] - b['socket_radius'])


def face_plane(doc):
    """Return (z0, slope) for the flat face z = z0 + slope*x, fit through the
    exit points of the two crown_ref holes. crown = perpendicular (-z) depth
    from the bore wall at that hole."""
    cr = doc['face']['crown_ref']
    pts = []
    for hid, crown in cr.items():
        x = hole(doc, hid)['position']
        pts.append((x, -(r_bore(doc, x) + crown)))
    (x1, z1), (x2, z2) = pts[0], pts[1]
    slope = (z2 - z1) / (x2 - x1)
    return z1 - slope * x1, slope


def face_z(doc, x):
    z0, slope = face_plane(doc)
    return z0 + slope * x


def hole_exit(doc, hid):
    """Geometry of a tonehole's angled chimney. Returns a dict:
      center  [x,y,z] where the chimney axis meets the face plane
      normal  [nx,ny,nz] outward chimney axis direction (unit)
      radius  hole radius (mm)
      tilt_deg
      path    axial chimney length, bore wall -> face (mm)
    Solved from the tilt and the face plane; nothing stored redundantly."""
    h = hole(doc, hid)
    x = h['position']
    a = h['diameter'] / 2.0
    th = math.radians(h['tilt_deg'])
    z0, slope = face_plane(doc)
    # chimney axis from bore centre (x,0,0) along d=(sin th, 0, -cos th);
    # intersect z = z0 + slope*x  ->  solve for parameter t (axis length)
    denom = math.cos(th) + slope * math.sin(th)
    t_axis = -(z0 + slope * x) / denom               # centre-line to plane
    cx = x + t_axis * math.sin(th)
    cz = -t_axis * math.cos(th)
    return {'center': [cx, 0.0, cz],
            'normal': [math.sin(th), 0.0, -math.cos(th)],
            'radius': a, 'tilt_deg': h['tilt_deg'],
            'path': t_axis - r_bore(doc, x)}


def manifest(doc):
    """Named openings for FEM face-tagging and CAD, all in one place.
    Coordinates are exact; the FEM tags an imported STEP face by nearest
    center match, so this is the contract between CAD and FEM."""
    z0, slope = face_plane(doc)
    return {
        'axes': 'x along bore from window (x=0); bore axis z=0; face z<0',
        'window': {'plane_x': 0.0, 'radius': doc['bore']['socket_radius']},
        'foot': {'center': [doc['bore']['length'], 0.0, 0.0],
                 'normal': [1.0, 0.0, 0.0], 'radius': doc['bore']['foot_radius']},
        'face': {'z0': z0, 'slope': slope},
        'holes': {hid: hole_exit(doc, hid) for hid in HOLE_IDS},
    }


# ---------------------------------------------------------- within-hand gaps
def within_hand_gaps(doc):
    """Exterior centre-to-centre distances for the constrained hole pairs,
    measured between chimney EXIT centres on the face (what the fingers feel)."""
    ex = {hid: hole_exit(doc, hid)['center'] for hid in HOLE_IDS}
    out = {}
    for a, b in doc['constraint']['within_hand_pairs']:
        dx = ex[b][0] - ex[a][0]
        dz = ex[b][2] - ex[a][2]
        out['%s-%s' % (a, b)] = math.hypot(dx, dz)
    return out


# ----------------------------------------------------------------- openwind
def to_openwind(doc, mask):
    """1-D transmission-line inputs for a fingering mask {hid: 'open'|'closed'}.
    Returns plain data (main_bore, holes_table, chart, temperature_c); the
    caller builds the ImpedanceComputation. Chimney length = the true angled
    path from hole_exit; open holes get the fitted inner correction dt_open."""
    b = doc['bore']
    main_bore = [[0.0, b['socket_radius'] / 1000.0],
                 [b['body_start'] / 1000.0, b['socket_radius'] / 1000.0],
                 [b['length'] / 1000.0, b['foot_radius'] / 1000.0]]
    dt = doc.get('hole_correction_dt_open', 0.0)
    rows = [['label', 'position', 'radius', 'chimney']]
    for hid in HOLE_IDS:
        h = hole(doc, hid)
        chim = hole_exit(doc, hid)['path'] + (dt if mask[hid] == 'open' else 0.0)
        rows.append([hid, h['position'] / 1000.0, h['diameter'] / 2000.0,
                     chim / 1000.0])
    labels = [ow_label(n) for n in doc['fingerings']]
    chart = [['label'] + labels]
    for hid in HOLE_IDS:
        chart.append([hid] + [_STATE[doc['fingerings'][n][hid]]
                              for n in doc['fingerings']])
    return main_bore, rows, chart, doc['target']['temperature_c']


def analytic_air_volume(doc, n=4000):
    """Air volume (mm^3) by integrating the bore + adding chimney frusta.
    The regen invariant compares this against the meshed/CAD solid so a broken
    generator fails loudly instead of shipping plausible-wrong geometry."""
    b = doc['bore']
    L = b['length']
    dx = L / n
    vol = 0.0
    for i in range(n):
        x = (i + 0.5) * dx
        vol += math.pi * r_bore(doc, x) ** 2 * dx
    for hid in HOLE_IDS:
        e = hole_exit(doc, hid)
        vol += math.pi * e['radius'] ** 2 * e['path']
    return vol


# ----------------------------------------------------------------- validate
def validate(doc):
    errs = []
    b = doc['bore']
    if not (b['body_start'] < b['length']):
        errs.append('bore.body_start must be < bore.length')
    ids = [h['id'] for h in doc['holes']]
    if ids != HOLE_IDS:
        errs.append('holes must be exactly %s in order, got %s'
                    % (HOLE_IDS, ids))
    for h in doc['holes']:
        if not (b['body_start'] < h['position'] < b['length']):
            errs.append('%s position %.1f outside body bore (%.1f..%.1f)'
                        % (h['id'], h['position'], b['body_start'], b['length']))
        if not (2.0 <= h['diameter'] <= 16.0):
            errs.append('%s diameter %.1f implausible' % (h['id'], h['diameter']))
        if abs(h['tilt_deg']) >= 75.0:
            errs.append('%s tilt %.1f too steep (>=75 deg)'
                        % (h['id'], h['tilt_deg']))
    for hid in doc['face']['crown_ref']:
        if hid not in ids:
            errs.append('face.crown_ref references unknown hole %s' % hid)
    for n in doc['fingerings']:
        m = doc['fingerings'][n]
        if sorted(m) != sorted(HOLE_IDS):
            errs.append('fingering %s does not cover all holes' % n)
        for st in m.values():
            if st not in _STATE:
                errs.append('fingering %s bad state %r' % (n, st))
    # chimney paths must be positive (face must clear the bore wall everywhere)
    for hid in HOLE_IDS:
        if hole_exit(doc, hid)['path'] <= 0.5:
            errs.append('%s chimney path <=0.5mm (face cuts into bore)' % hid)
    if errs:
        raise ValueError('spec invalid:\n  ' + '\n  '.join(errs))
    return doc


def with_holes(doc, positions=None, diameters=None, tilts=None, length=None):
    """Return a copy with solver overrides applied (for the re-solve loop)."""
    d = copy.deepcopy(doc)
    if length is not None:
        d['bore']['length'] = length
    for h in d['holes']:
        if positions and h['id'] in positions:
            h['position'] = positions[h['id']]
        if diameters and h['id'] in diameters:
            h['diameter'] = diameters[h['id']]
        if tilts and h['id'] in tilts:
            h['tilt_deg'] = tilts[h['id']]
    return d
