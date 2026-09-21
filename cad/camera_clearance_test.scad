// Dovetail clearance coupons, to pick dt_cl for a given printer before
// printing the carriage. 0 prints the three side by side; 1, 2, 3 prints one
// (0.20, 0.35 and 0.50 mm).
include <camera_mount.scad>

variant = 0;

if (variant == 0) base_jeu_all();
else              base_jeu_labeled(dt_cl_test_values[variant - 1]);
