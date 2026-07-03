# Build the M2 folded low D body (foot U-fold) in a build123d-mcp session.
#
# Run AFTER loading foldflute/prelude.py into the session. The parameter
# defaults below are the solved values from reports/solved_designs.json
# (m2_folded); regenerate them with optimize_sweep()/solve_foot_fold() if the
# spec changes. All dimensions mm; body coordinates: x from the body start
# (headjoint tube end, arc 120), fold in the XY plane, hole exits +z.
#
# Print architecture (per reports/REPORT.md):
#   section_a     x 0..172 + 20 tenon; monolithic, printed vertically
#                 (keeps the sealing spigot and tenon round); carries h1.
#   left_shell    front (+z) half of section B; carries h2..h6 entirely.
#   right_shell   back (-z) half of section B.
#   alignment_pin loose dowel (print 3 + spares).
from build123d import *
import math

# ---- solved acoustics (window-arc mm) -> reports/solved_designs.json ----
L_END   = 574.567          # bore end from window (D4 solve, bend-corrected)
S_BEND  = 464.895          # bend start from window
BODY0   = 120.0            # window -> body socket (headjoint tube end)
BEND_R  = 33.0             # 1.5 x bore diameter (legal minimum)
POSITIONS = {'h1': 283.62, 'h2': 319.27, 'h3': 358.91,
             'h4': 411.12, 'h5': 431.90, 'h6': 488.99}
DIAMETERS = {'h1': 11.0, 'h2': 10.5, 'h3': 10.0,
             'h4': 10.0, 'h5': 10.5, 'h6': 11.0}

# ---- body parameters ----
BORE_R   = 11.0
OD_R     = 14.3            # body OD 28.6 = sleeve OD
SPIG_R   = 12.4            # spigot OD 24.8 into the 25.0 ID sleeve
SPIG_LEN = 56.0
JOINT_X  = 172.0           # A/B joint (between h1 and h2), body x
TENON_R, TENON_LEN = 13.0, 20.0
SOCK_R, SOCK_DEPTH = 13.1, 20.4
RING_R, RING_X0, RING_X1 = 15.4, 172.0, 194.0
BOSS_R, BOSS_H = 4.5, 5.2
BOSS_XY = [(215.0, 16.0), (300.0, -16.0), (394.2, 33.0)]
PIN_HOLE_R, PIN_HOLE_DEPTH = 1.25, 2.2
BIG = 1000.0

SBEND = S_BEND - BODY0
TAIL = L_END - S_BEND - math.pi * BEND_R
FOOT_X = SBEND - TAIL
PASS_Y = 2.0 * BEND_R

hole_geo = {}
for hid, s_arc in POSITIONS.items():
    s = s_arc - BODY0
    if s <= SBEND:
        hx, hy = s, 0.0
    else:
        ph = (s - SBEND) / BEND_R
        hx, hy = SBEND + BEND_R * math.sin(ph), BEND_R * (1.0 - math.cos(ph))
    hole_geo[hid] = (hx, hy, DIAMETERS[hid] / 2.0)

# ---- canonical air column ----
with BuildLine() as bore_path:
    Line((0, 0, 0), (SBEND, 0, 0))
    JernArc(start=(SBEND, 0), tangent=(1, 0), radius=BEND_R, arc_size=180)
    Line((SBEND, PASS_Y, 0), (FOOT_X, PASS_Y, 0))
with BuildSketch(Plane.YZ) as bore_prof:
    Circle(BORE_R)
air = sweep(bore_prof.sketch, path=bore_path.line)
show(air, 'air_solid')

# ---- canonical body: spigot + swept tube + joint ring - air - holes ----
with BuildLine() as body_path:
    Line((SPIG_LEN, 0, 0), (SBEND, 0, 0))
    JernArc(start=(SBEND, 0), tangent=(1, 0), radius=BEND_R, arc_size=180)
    Line((SBEND, PASS_Y, 0), (FOOT_X, PASS_Y, 0))
with BuildSketch(Plane.YZ.offset(SPIG_LEN)) as body_prof:
    Circle(OD_R)
main_tube = sweep(body_prof.sketch, path=body_path.line)
spigot = Cylinder(SPIG_R, SPIG_LEN, rotation=(0, 90, 0)).move(
    Location((SPIG_LEN / 2.0, 0, 0)))
ring = Cylinder(RING_R, RING_X1 - RING_X0, rotation=(0, 90, 0)).move(
    Location(((RING_X0 + RING_X1) / 2.0, 0, 0)))
body = spigot.fuse(main_tube).fuse(ring) - air
cutters = None
for hid in ['h1', 'h2', 'h3', 'h4', 'h5', 'h6']:
    hx, hy, hr = hole_geo[hid]
    c = Cylinder(hr, 14.0).move(Location((hx, hy, 13.0)))
    cutters = c if cutters is None else cutters.fuse(c)
body = body - cutters
for bx, by in BOSS_XY:                     # seam alignment bosses
    body = body.fuse(Cylinder(BOSS_R, BOSS_H).move(Location((bx, by, 0))))
show(body, 'body_solid')

# ---- section A (monolithic) and section B (shell pair) ----
cutA = Box(BIG, BIG, BIG, align=(Align.MAX, Align.CENTER, Align.CENTER)).move(
    Location((JOINT_X, 0, 0)))
cutB = Box(BIG, BIG, BIG, align=(Align.MIN, Align.CENTER, Align.CENTER)).move(
    Location((JOINT_X, 0, 0)))
tenon = (Cylinder(TENON_R, TENON_LEN, rotation=(0, 90, 0))
         - Cylinder(BORE_R, TENON_LEN, rotation=(0, 90, 0))).move(
    Location((JOINT_X + TENON_LEN / 2.0, 0, 0)))
socket_recess = Cylinder(SOCK_R, SOCK_DEPTH, rotation=(0, 90, 0)).move(
    Location((JOINT_X + SOCK_DEPTH / 2.0, 0, 0)))
sectionA = (body & cutA).fuse(tenon)
sectionB = (body & cutB) - socket_recess
show(sectionA, 'section_a')
show(sectionB, 'section_b')

halfF = Box(BIG, BIG, BIG, align=(Align.CENTER, Align.CENTER, Align.MIN))
halfK = Box(BIG, BIG, BIG, align=(Align.CENTER, Align.CENTER, Align.MAX))
shell_front = sectionB & halfF
shell_back = sectionB & halfK
for bx, by in BOSS_XY:                     # dowel recesses in both seam faces
    shell_front = shell_front - Cylinder(PIN_HOLE_R, PIN_HOLE_DEPTH).move(
        Location((bx, by, PIN_HOLE_DEPTH / 2.0)))
    shell_back = shell_back - Cylinder(PIN_HOLE_R, PIN_HOLE_DEPTH).move(
        Location((bx, by, -PIN_HOLE_DEPTH / 2.0)))
show(shell_front, 'left_shell')
show(shell_back, 'right_shell')
show(Cylinder(1.15, 4.0), 'alignment_pin')

# ---- print orientations (z = bed) ----
show(sectionA.rotate(Axis.Y, 90).move(Location((0, 0, JOINT_X + TENON_LEN))),
     'section_a_print')                        # spigot up, tenon on bed
show(shell_back.rotate(Axis.X, 180), 'right_shell_print')  # seam face down
# left_shell already sits seam-face-down.
