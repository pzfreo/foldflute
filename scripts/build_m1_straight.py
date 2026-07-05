# Build the M1 straight low D body (two vertical-print barrel sections) in a
# build123d-mcp session. Run AFTER loading foldflute/prelude.py. Parameter
# defaults are the solved values from reports/solved_designs.json
# (m1_straight). All dimensions mm; x from the body start (arc 120).
from build123d import *

# ---- solved acoustics (window-arc mm) ----
L_END = 571.603
BODY0 = 120.0
POSITIONS = {'h1': 283.62, 'h2': 319.27, 'h3': 358.92,
             'h4': 411.13, 'h5': 431.89, 'h6': 488.34}
DIAMETERS = {'h1': 11.0, 'h2': 10.5, 'h3': 10.0,
             'h4': 10.0, 'h5': 10.5, 'h6': 11.0}

# ---- body parameters (shared with the M2 design) ----
BORE_R, OD_R = 11.0, 14.3
SPIG_R, SPIG_LEN = 12.4, 56.0
JOINT_X = 219.0            # between h2 (edge 204.5) and h3 (edge 233.9)
TENON_R, TENON_LEN = 13.0, 20.0
SOCK_R, SOCK_DEPTH = 13.1, 20.4
RING_R, RING_X0, RING_X1 = 15.4, 219.0, 241.0
BIG = 1000.0
ARC = L_END - BODY0

spigot = Cylinder(SPIG_R, SPIG_LEN, rotation=(0, 90, 0)).move(
    Location((SPIG_LEN / 2.0, 0, 0)))
tube = Cylinder(OD_R, ARC - SPIG_LEN, rotation=(0, 90, 0)).move(
    Location(((SPIG_LEN + ARC) / 2.0, 0, 0)))
ring = Cylinder(RING_R, RING_X1 - RING_X0, rotation=(0, 90, 0)).move(
    Location(((RING_X0 + RING_X1) / 2.0, 0, 0)))
air = Cylinder(BORE_R, ARC, rotation=(0, 90, 0)).move(
    Location((ARC / 2.0, 0, 0)))
body = spigot.fuse(tube).fuse(ring) - air
cutters = None
for hid in ['h1', 'h2', 'h3', 'h4', 'h5', 'h6']:
    c = Cylinder(DIAMETERS[hid] / 2.0, 14.0).move(
        Location((POSITIONS[hid] - BODY0, 0, 13.0)))
    cutters = c if cutters is None else cutters.fuse(c)
body = body - cutters
show(body, 'm1_body')

cutA = Box(BIG, BIG, BIG, align=(Align.MAX, Align.CENTER, Align.CENTER)).move(
    Location((JOINT_X, 0, 0)))
cutB = Box(BIG, BIG, BIG, align=(Align.MIN, Align.CENTER, Align.CENTER)).move(
    Location((JOINT_X, 0, 0)))
tenon = (Cylinder(TENON_R, TENON_LEN, rotation=(0, 90, 0))
         - Cylinder(BORE_R, TENON_LEN, rotation=(0, 90, 0))).move(
    Location((JOINT_X + TENON_LEN / 2.0, 0, 0)))
socket_recess = Cylinder(SOCK_R, SOCK_DEPTH, rotation=(0, 90, 0)).move(
    Location((JOINT_X + SOCK_DEPTH / 2.0, 0, 0)))
m1_a = (body & cutA).fuse(tenon)
m1_b = (body & cutB) - socket_recess
show(m1_a, 'm1_section_a')
show(m1_b, 'm1_section_b')

# ---- print orientations (both vertical; brim recommended) ----
show(m1_a.rotate(Axis.Y, 90).move(Location((0, 0, JOINT_X + TENON_LEN))),
     'm1_a_print')                             # spigot up, tenon on bed
show(m1_b.rotate(Axis.Y, -90).move(Location((0, 0, -JOINT_X))),
     'm1_b_print')                             # socket ring on bed, foot up
