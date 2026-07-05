"""foldflute solid generation — build123d, the ONE and only CAD path.

Runs in .venv-cad (build123d). Every solid derives from foldflute.geometry
reading the canonical doc, so the CAD and the 1-D/FEM models cannot describe
different geometry. `export_air` writes the exact solid the FEM gate meshes,
plus manifest.json (from geometry.manifest) that names its openings — that
manifest is the contract that lets the gate run on the literal exported STEP
instead of a reconstruction.

Convention matches geometry.py: x along the bore (window at 0), bore axis on
z=0, flat front face below at z<0.
"""
import json
import math
import os

from build123d import (Align, Axis, Box, Cone, Cylinder, Location,
                       export_step)

from . import geometry as G

_BIG = 4000.0


def _taper_bore(doc, radius_fn, extra_foot=0.0):
    """Air/solid of revolution following the bore: cylinder to body_start then
    conical frusta to the foot. radius_fn(x) supplies the radius."""
    b = doc['bore']
    parts = [Cylinder(radius_fn(0.0), b['body_start'], rotation=(0, 90, 0))
             .move(Location((b['body_start'] / 2.0, 0, 0)))]
    nseg = 12
    x0 = b['body_start']
    span = b['length'] + extra_foot - x0
    for i in range(nseg):
        xa = x0 + span * i / nseg
        xb = x0 + span * (i + 1) / nseg
        ra, rb = radius_fn(xa), radius_fn(xb)
        if abs(ra - rb) < 1e-6:
            seg = Cylinder(ra, xb - xa, rotation=(0, 90, 0))
        else:
            seg = Cone(ra, rb, xb - xa, rotation=(0, 90, 0))
        parts.append(seg.move(Location(((xa + xb) / 2.0, 0, 0))))
    out = parts[0]
    for p in parts[1:]:
        out = out + p
    return out


def _face_halfspace(doc, keep='below', pad=0.0):
    """A big box whose top face lies on the tilted face plane (offset by pad
    along +z). keep='below' returns material on the bore side (z<face)."""
    z0, slope = G.face_plane(doc)
    ang = math.degrees(math.atan(slope))
    align = (Align.CENTER, Align.CENTER,
             Align.MAX if keep == 'below' else Align.MIN)
    box = Box(_BIG, _BIG, _BIG, align=align)
    box = box.rotate(Axis.Y, ang)
    return box.move(Location((0, 0, z0 + pad)))


def build_air(doc):
    """The air column: tapered bore + angled tonehole chimneys, cut flush at
    the face plane so each chimney ends in a coplanar exit disc."""
    air = _taper_bore(doc, lambda x: G.r_bore(doc, x))
    for hid in G.HOLE_IDS:
        e = G.hole_exit(doc, hid)
        x = G.hole(doc, hid)['position']
        th = e['tilt_deg']
        length = e['path'] + G.r_bore(doc, x) + 6.0     # from axis, past face
        chim = (Cylinder(e['radius'], length, align=(Align.CENTER, Align.CENTER,
                                                     Align.MIN))
                .rotate(Axis.Y, 180.0 - th)             # point outward (-z side)
                .move(Location((x, 0, 0))))
        air = air + chim
    air = air & _face_halfspace(doc, keep='below')      # trim to the face
    return air


def build_body(doc, wall=2.5, back_depth=None):
    """A minimal printable-shaped body: tapered outer wall following the bore
    +wall, flat front face, bored out by the air. Not the final print part
    (no joint/sections/chamfers) — enough to view proportions and to run the
    as-exported gate on a real closed solid."""
    b = doc['bore']
    if back_depth is None:
        back_depth = b['socket_radius'] + wall + 12.0
    outer = _taper_bore(doc, lambda x: G.r_bore(doc, x) + wall, extra_foot=0.0)
    # add a rectangular back so the body is a solid slab behind the bore
    z0, slope = G.face_plane(doc)
    slab = (Box(b['length'] - b['body_start'], 2 * (b['socket_radius'] + wall),
                back_depth, align=(Align.MIN, Align.CENTER, Align.MAX))
            .move(Location((b['body_start'], 0, b['socket_radius'] + wall))))
    body = outer + slab
    body = body & _face_halfspace(doc, keep='below', pad=0.0)
    body = body - build_air(doc)
    return body


def export_air(doc, out_dir):
    """Write air.step (the FEM gate's input) + manifest.json (its face tags)."""
    os.makedirs(out_dir, exist_ok=True)
    air = build_air(doc)
    step = os.path.join(out_dir, 'air.step')
    export_step(air, step)
    man = G.manifest(doc)
    man['air_volume_mm3'] = air.volume
    with open(os.path.join(out_dir, 'manifest.json'), 'w') as fh:
        json.dump(man, fh, indent=2)
    return step, air.volume


def export_body(doc, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    body = build_body(doc)
    step = os.path.join(out_dir, 'body.step')
    export_step(body, step)
    return step, body.volume
