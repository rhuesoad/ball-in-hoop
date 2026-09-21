// Whole bench, for checking fits. Not a print file.
//   "assembled", "exploded"          : the hoops, motor support and link
//   "camera"                         : rail and carriage on the M6 grid
//   "drawing_motor" / "_outer" / "_inner" : three-view 2D drawings
include <bench.scad>
include <camera_mount.scad>

view    = "assembled";
explode = 40;
car_y   = 0;

if      (view == "assembled")     orienter_axes() assembly(explode = 0);
else if (view == "exploded")      orienter_axes() assembly(explode = explode);
else if (view == "drawing_motor") mise_en_plan_3_vues() motor_support_with_fixations();
else if (view == "drawing_outer") mise_en_plan_3_vues() outer_hoop();
else if (view == "drawing_inner") mise_en_plan_3_vues() inner_hoop();
else if (view == "camera")
    rotate([0, 0, 180]) translate([0, 300, -height_hole]) {
        color("Tan")       rail();
        color("Steelblue") carriage(car_y);
        grid_ref();
    }
