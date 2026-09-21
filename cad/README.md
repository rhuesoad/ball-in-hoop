# CAD — printed parts

OpenSCAD source of every printed part (tested on OpenSCAD 2021.01). Open a part
file, render (F6), export to STL.

| File | Part |
|------|------|
| `outer_hoop_down.scad`, `outer_hoop_up.scad` | Outer hoop, two halves |
| `inner_hoop_down.scad`, `inner_hoop_up.scad` | Inner hoop, two halves |
| `motor_support.scad` | Motor support with its floor brackets |
| `motor_link.scad` | Link between the motor shaft and the hoops |
| `odrive_plate.scad` | Plate holding the ODrive S1 |
| `camera_rail.scad` | Camera rail, bolted to the M6 grid |
| `camera_carriage.scad` | Camera carriage, slides on the rail |
| `camera_clearance_test.scad` | Coupons to tune the rail clearance to a printer |
| `assembly.scad` | Whole bench, exploded view and 2D drawings (not a print file) |

The part files only select what to render. The geometry lives in `bench.scad`
(hoops and motor mount) and `camera_mount.scad` (camera rail).

`lib/nutsnbolts.scad` is "Norm Nuts and Bolts" by Johannes Kneer, under GPLv3.
