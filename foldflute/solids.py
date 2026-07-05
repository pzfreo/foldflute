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
                       Plane, Vector, export_step)

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


def _face_halfspace(doc, keep='bore', pad=0.0):
    """Half-space solid bounded exactly by the tilted face plane z = z0+slope*x.
    keep='bore' returns the material side (toward the bore, +normal); 'outside'
    the room side. Built in a Plane frame whose Z axis is the face normal, so
    there is no rotation-order ambiguity."""
    z0, slope = G.face_plane(doc)
    # bore-side normal points away from the face into the material (+z-ish)
    n = Vector(-slope, 0.0, 1.0)
    n = n / n.length
    pl = Plane(origin=Vector(0, 0, z0) + n * pad, z_dir=n)
    zalign = Align.MIN if keep == 'bore' else Align.MAX
    box = Box(_BIG, _BIG, _BIG, align=(Align.CENTER, Align.CENTER, zalign))
    return pl * box


def _face_start(doc):
    """x where the flat face begins. Above this (toward the fipple) the body
    is a full round tube — a single tilted plane through the hole exits would
    otherwise rise above the bore bottom and slice the head bore."""
    return G.hole(doc, 'h1')['position'] - 14.0


def _keep_region(doc, pad=0.0):
    """Material region = full head tube (x < face_start) UNION the face
    half-space over the hole span. Intersecting with this trims chimneys at
    the face without eating the round head bore."""
    xs = _face_start(doc)
    head = Box(2 * xs + _BIG, _BIG, _BIG,
               align=(Align.MAX, Align.CENTER, Align.CENTER)).move(
        Location((xs, 0, 0)))
    return head + _face_halfspace(doc, keep='bore', pad=pad)


def build_air(doc):
    """The air column: tapered bore + angled tonehole chimneys, cut flush at
    the face plane over the hole span (round head bore preserved)."""
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
    air = air & _keep_region(doc)                       # trim chimneys to face
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
    body = body & _face_halfspace(doc, keep='bore', pad=0.0)
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
