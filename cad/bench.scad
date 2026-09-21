include <lib/nutsnbolts.scad>

// Set once in the original single file, and so applied to every part in it.
$fn = 48;

//----------------------------------------------------------------------------------
// ------------ Parameters ------------ 
//----------------------------------------------------
R = 1.5; // radius of the o-ring core
r2 = (204+2*R)/2; // inner radius of the hoop (from center of the hoop to the center of the o ring profile)
r3 = 84/2+R+1.5; // outter radius of the inner hoop 
h1 = 4; // height of the base of the hoop
h2 = 8; // distance from the base to the center of the o-ring profile
h3 = 14.7; // distance between the o-rings (center to center)
h4 = 2.5; //distance from the outter o-ring to the outter edge of the hoop

ushape_oRings_dist = 11; // distance between the center of the inner and outer o-ring on the ushape inner hoop
ushape_R = 50;

h_tot = h1+h2+h3+h4;
R_ext = 0.25;
r_r = 2; // depth of the rim
r1 = r2+2*R+2+10; // outter radius of the hoop

echo(h_tot+h1);

r_in = 27.5;


echo("External diameter: ", 2*r1);
echo("Diam to outer o-ring: ", 2*r2-3);
echo("Diam to outer inner o-ring: ", 2*(ushape_R+ushape_oRings_dist/2));
echo("Diam to inner inner o-ring: ", 2*(ushape_R-ushape_oRings_dist/2));
echo("Inner diameter of the outer hoop is: ", 2*(r2+R-dy));

//echo("Outter diameter of the hoop is: ", 2*r1);
//echo("Outter diameter of the inner hoop is: ", 2*r3);
//echo("Inner diameter of the outer hoop is: ", 2*(r2+R-dy));
//echo("Diameter of the inner part of the inner ushape (to the center of the o-ring)", 2*(ushape_R-ushape_oRings_dist/2));
//echo("Diameter of the outer part of the inner ushape (to the center of the o-ring)", 2*(ushape_R+ushape_oRings_dist/2));

d_nevim = 0.5;
h_nevim = 0.5;

d_rib = 16;

//d_tooth = 0.1231;
d_tooth = 0.25;

alpha = 55;
beta = 25;
gamma = 50;

// Calculate some values for a bit (not much actually) more clear calculations in the following script 
hx = h1+h2+R*sin(alpha);
dx = r1-r2-R*(1+cos(alpha));
hy = h2-(R+R_ext)*sin(beta);
dy = (R+R_ext)*cos(beta);

hz = r_r*tan(gamma);

r_groove = r1-dx-sqrt(d_nevim*d_nevim/2)+r_r; //thickness of the wall between the outter side of the large hoop and the bottom of the groove

part1_height = hx;
part2_height = 2*(h3/2+h2+h1-hx);
part3_height = h1+h2+h3+h4 - part1_height - part2_height;


//----------------------------------------------------
// ------------ Modelling ------------ 
//----------------------------------------------------
// outter o-rings
module outter_oring_lower() {
    translate([r2+R,h1+h2,0]){
        circle(R, $fn=50);
    };
};
module outter_oring_upper() {
    translate([0,h3,0]){
        outter_oring_lower();
    };
};

// Mounting holes
module mounting_hole_m3_tall(){
    cylinder(h=h1+h2+h3+h4+1, r=1.6, center=false,$fn=50);
    cylinder(h=1.5, r=5.5/2, center=false, $fn=50);
}
module mounting_hole_m3(){
    cylinder(h=h1+h2+h3+h4+1, r=1.6, center=false,$fn=50);
}
module motor_hole(){
    // Passage M4 (Ø4 mm + jeu 0.4 mm)
    cylinder(h=h1, r=2.2, center=false, $fn=50);
    // Lamage tête M4 (Ø7 mm, hauteur tête = 4 mm)
    cylinder(h=0.5, r=3.5, center=false, $fn=50);
}

// ------------ outter hoop (the larger one) ------------ 

//// Profiles of the parts of the outter hoop
// Profile 1 - the bottom one
//      before the subtraction of the o-rings
module profile_part01() {
    polygon(points=[[r2+R-dy,0],[r1,0],[r1,hx],[r1-dx,hx], [r2+R,h1+h2],[r2+R-dy,h1+hy],[r2+R-dy,h1]]);
};

// model of the tooth biting into the oring (parameter d_tooth)
module profile_part01_tooth() {
    polygon(points=[[r2+R-R*cos(beta),h1+h2-R*sin(beta)],[r2+R-R*cos(beta)+d_tooth,h1+h2-R*sin(beta)],[r2+R-R*cos(beta)+d_tooth,h1+h2-R],[r2+R-R*cos(beta),h1+h2-R]]);
}

//      after the subtraction of the o-ring
module profile_part1() {
    union() {
        difference() {
            profile_part01();
            outter_oring_lower();
        };
        profile_part01_tooth();
    }
};

// Profile 2 - the middle one
module profile_part2_half() {
    polygon(points=[[r1,hx], [r1-dx,hx], [r1-dx-sqrt(d_nevim*d_nevim/2),hx+sqrt(d_nevim*d_nevim/2)],[r1-dx-sqrt(d_nevim*d_nevim/2),hx+sqrt(d_nevim*d_nevim/2)+h_nevim], [r_groove, hx+sqrt(d_nevim*d_nevim/2)+h_nevim+hz],[r_groove, h1+h2+h3/2],[r1,h1+h2+h3/2]]);
};
module profile_part2() {
    union(){
        profile_part2_half();
        translate([0,2*(h1+h2+h3/2),0]) {
            mirror([0,1,0]){
                color("Cyan") profile_part2_half();
            };
        };
    };
};

// Profile 3 - the bottom one
//      before the subtraction of the o-rings
module profile_part03() {
    polygon(points=[[r1,h1+h2+h3-R*sin(alpha)],[r1-dx,h1+h2+h3-R*sin(alpha)], [r2+R,h1+h2+h3], [r2+R-dy,h1+h2+h3+(R+R_ext)*sin(beta)], [r2+R-dy,h1+h2+h3+h4], [r1,h1+h2+h3+h4]]);
};

// model of the tooth biting into the oring (parameter d_tooth)
module profile_part03_tooth() {
    polygon(points=[[r2+R-R*cos(beta),h1+h2+h3+R*sin(beta)],[r2+R-R*cos(beta)+d_tooth,h1+h2+h3+R*sin(beta)],[r2+R-R*cos(beta)+d_tooth,h1+h2+h3+R],[r2+R-R*cos(beta),h1+h2+h3+R]]);
}

//      after the subtraction of the o-ring
module profile_part3() {
    union() {
        difference() {
            profile_part03();
            outter_oring_upper();
        }
        profile_part03_tooth();
    }
};


//// Go to 3D
// Outter mounting holes
module outter_mounting_holes() {
    for (i=[1:12]) {
        rotate([0,0,i*30])
            translate([0,r1-(r1-r_groove)/2+0.3,0])
                mounting_hole_m3_tall();
    };
};

// Revolve the profiles and subtract the outter mounting holes
module part1(){
    union() {
        difference(){
            rotate_extrude($fn = 360){
                profile_part1();
            };
            outter_mounting_holes();
        };
        base();
    };
};

module part2(){
    difference(){
        rotate_extrude($fn = 360){
            profile_part2();
        };
        outter_mounting_holes();
    }
};

module part3(){
    difference(){
        rotate_extrude($fn = 360){
            profile_part3();
        };
        outter_mounting_holes();
    }
};

//// Test parts
module part1_test(revolve_angle){
    union() {
        difference(){
            rotate(90-revolve_angle/2) {
                rotate_extrude(angle=revolve_angle, $fn = 360){
                    profile_part1();
                };
            }
            outter_mounting_holes();
        };
    };
};

module part2_test(revolve_angle){
    difference(){
        rotate(90-revolve_angle/2) {
            rotate_extrude(angle=revolve_angle, $fn = 360){
                profile_part2();
            };
        }
        outter_mounting_holes();
    }
};

module part3_test(revolve_angle){
    difference(){
        rotate(90-revolve_angle/2) {
            rotate_extrude(angle=revolve_angle, $fn = 360){
                profile_part3();
            };
        }
        outter_mounting_holes();
    }
};

// ------------ The base ------------ 
// inner mounting holes
module inner_mounting_holes(){
    for (i=[1:6]) {
        rotate([0,0,i*60])
            translate([0,ushape_R,0]) {
                mounting_hole_m3();
                //scale(1.05)
                //    translate([0,0,2])
                //        nut("M3");
                // Lamage cylindrique Ø5 mm pour noyer la tête de vis M3 (hauteur tête = 1.5 mm)
                cylinder(h=1.5, r=5.5/2, center=false, $fn=50);
            };
    }

    // Bldc mounting holes - old gimbal motor
    //for (i=[1:4]) {
    //    rotate([0,0,i*90])
    //        translate([0, 25/2,0])
    //            mounting_hole_m3();
    //    rotate([0,0,i*90+45])
    //        translate([0, 30/2,0])
    //            mounting_hole_m3();
    //}
    
    // Bldc mounting holes - new ODrive motor
    //for (i=[1:4]) {
    //    rotate([0,0,i*90+45])
    //        translate([0, 20/2,0])
    //            mounting_hole_m3();
    //}
    
    for (i=[1:4]){
        rotate([0,0,i*90])
            translate([16.7, 16.7 , 0])
                motor_hole();
    }

    cylinder(h=h1+1, r=7.2, center=false,$fn=50);
}


module inner_circle_small(height) {
    cylinder(h=height, r=r_in, center=false,$fn=180);
}

module inner_circle_large(height) {
    cylinder(h=height, r=r2+R-dy+1, center=false,$fn=180);
}


// Ribs joining the outter hoop to the inner circle
// Model one rib
module bottom_planar_rib() {
    translate([r_in-d_rib,0,0]){
        square([d_rib,(r2+R-dy)],false);
    };
};

// copy and rotate it four times to get half of the total ribs
module bottom_planar_4ribs() {
    for (i=[1:4]) {
        rotate([0,0,i*90])
            bottom_planar_rib();
    }
};

// mirror the four ribs to get the remaining four
module bottom_8ribs(){
    bottom_planar_4ribs();
    mirror([1,0,0]){
        bottom_planar_4ribs();
    };
};

// extrude all the eight ribs and take the intersetion of the extruded ribs with the larger inner circle (full base)
module ribs(h) {
    intersection(){
        linear_extrude(height = h) {
            bottom_8ribs();
        };
        inner_circle_large(h);
    };
};

module inner_mounting_circle() {
    linear_extrude(height = h1)
        difference() {
            circle(r=ushape_R+dy+ushape_oRings_dist/2, $fn=200);
            circle(r=ushape_R-dy-ushape_oRings_dist/2, $fn=200);
        }
}

// Model the base - take union of the ribs, inner circle and subtract the inner mounting holes
module base() {
    difference(){
        union(){
            ribs(h1);
            inner_circle_small(h1);
            inner_mounting_circle();
        }
        inner_mounting_holes();
    };
}

//----------------------------------------------------
// ------------------- Inner hoop  ------------------- 
//----------------------------------------------------
module profile_part01_inner() {
    translate([-r1+r2+r3,0,0])
        mirror([1,0,0])
            translate([-r1,0,0]) {
                difference(){
                    profile_part1();
                    polygon(points=[[r2+R-dy-0.1,-0.1],[r1+0.1,-0.1],[r1,h1],[r2+R-dy,h1]]);
                }
            }
};


module profile_part02_inner() {
    translate([-r1+r2+r3,0,0])
        mirror([1,0,0])
            translate([-r1,0,0]) {
                profile_part2();
            }
};

module profile_part03_inner() {
    translate([-r1+r2+r3,0,0])
        mirror([1,0,0])
            translate([-r1,0,0]) {
                profile_part3();
            }
};

//inner hoop mounting pillar
module inner_hoop_pillar(height, length=7, width=6.4, countersunk=false) {
    translate([0,47/2,0]) {
        difference() {
            union() {
                translate([0,0,h1])
                    cylinder(h=height, r=width/2, center=false, $fn=50);
                translate([-width/2,0,h1])
                    cube([width,width+length,height], center=false);
            }
            mounting_hole_m3();
            if(countersunk) {
                translate([0,0,height+h1-1.6])
                    cylinder(h=1.6,r1=1.6,r2=2.8, center=false, $fn=50);
            }        
        }

    }
}

// inner mounting holes
module inner_hoop_pillars(height, length=7, width=6.4, number=6, countersunk=false){
    for (i=[1:number]) {
        rotate([0,0,i*60])
            inner_hoop_pillar(height, length, width, countersunk);
    }
}	

module inner_hoop_part1() {
	union() {
		rotate_extrude($fn = 200)
	            profile_part01_inner();
	    inner_hoop_pillars(hx-h1);
	}
}
module inner_hoop_part2() {
	union() {
		rotate_extrude($fn = 200)
		            profile_part02_inner();
		translate([0,0,hx-h1])
			inner_hoop_pillars(part2_height);
		}
}
module inner_hoop_part3() {
	union() {
		rotate_extrude($fn = 200)
		            profile_part03_inner();
		translate([0,0,hx-h1+part2_height])
			inner_hoop_pillars(part3_height);
	}
}

module backBlackSheet() {
    difference() {
        circle(r=r1, center=true, $fn=200);
        
        // outer holes
        for (i=[1:12]) {
            rotate([0,0,i*30])
                translate([0,r1-(r1-r_groove)/2+0.3,0])
                    circle(r=1.5, center=true, $fn=200);
        }
        
        for (i=[1:6]) {
            rotate([0,0,i*60])
                translate([0,ushape_R,0]) {
                    circle(r=3, center=true, $fn=200);
                };
        }

//        // Bldc mounting holes
//        for (i=[1:4]) {
//            rotate([0,0,i*90])
//                translate([0, 25/2,0])
//                    circle(r=1.5, center=true, $fn=200);
//            rotate([0,0,i*90+45])
//                translate([0, 30/2,0])
//                    circle(r=1.5, center=true, $fn=200);
//        }
        
    // Bldc mounting holes - new ODrive motor
    for (i=[1:4]) {
        rotate([0,0,i*90+45])
            translate([0, 20/2,0])
                circle(r=1.5, center=true, $fn=200);
    }
    
    circle(r=7.5, center=false,$fn=50);
    
    }
}

//----------------------------------------------------
// --------------- U-shape inner hoop  ---------------
//----------------------------------------------------

module profile_outer() {
    profile_part1();
    profile_part2();
    profile_part3();
}

module half_profile_ushape() {
    translate([0,0,0])
    difference(){
        translate([-r2-R,0,0]){
            profile_outer();
        }
        translate([ushape_oRings_dist/2,0,0])
            square([30, 30]);
        translate([-30,0,0])
            square([60, h1]);
    }
}

module profile_ushape() {
    union() {
        translate([-ushape_oRings_dist/2,0,0])
            half_profile_ushape();
        mirror([1,0,0])
            translate([-ushape_oRings_dist/2,0,0])
                half_profile_ushape();
    }
}

//half_profile_ushape();
module ushape(rrr, angl=270) {
    difference() {
        rotate([0,0,135]) {
            // The main body
            rotate_extrude(angle=angl, $fn=200)
                translate([rrr, 0, 0])
                    profile_ushape();
            
            // Endings
            translate([rrr,0,0])
                rotate_extrude($fn = 200, angle=180) {
                    translate([-ushape_oRings_dist/2,0,0])
                        half_profile_ushape();
                }    
            
            rotate([0,0,angl])
                translate([rrr,0,0])
                    rotate([0,0,180])
                    rotate_extrude($fn = 200, angle=180) {
                        translate([-ushape_oRings_dist/2,0,0])
                            half_profile_ushape();
                    }
        }
        
        inner_mounting_holes();
    }
}


//            circle(r=ushape_R+r1-r2+R-ushape_oRings_dist/2, $fn=200);


// ============================================================================
// DÉCOUPE HORIZONTALE — plan de coupe à mi-hauteur entre les deux O-rings
// z_cut = h1 + h2 + h3/2
// ============================================================================
z_cut = h1 + h2 + h3/3;

// Boîte de découpe générique — assez grande pour englober n'importe quelle pièce
bbox = 2 * r1 ;

module outer_hoop() {
    union(){part1();part2();part3();}
}

module outer_down() {
    intersection() {
        outer_hoop();
        translate([-bbox, -bbox, 0]){
            cube([bbox*2, bbox*2, z_cut]);
        }
        
    }
}

module outer_up() {
    intersection() {
        outer_hoop();
        translate([-bbox, -bbox, z_cut]){
            cube([bbox*2, bbox*2, z_cut]);
        }
        
    }
}



//translate([0, 0, 2*z_cut]) rotate([0, 180, 0]) color("Cyan") outer_up();

module inner_hoop() {
    ushape(ushape_R, 270);
}

module inner_down() {
    intersection() {
        inner_hoop();
        translate([-bbox, -bbox, 0]){
            cube([bbox*2, bbox*2, z_cut]);
        }
    }
}

module inner_up() {
    intersection() {
        inner_hoop();
        translate([-bbox, -bbox, z_cut]){
            cube([bbox*2, bbox*2, z_cut]);
        }
    }
}

//color("Pink") translate([300, 0, -h1]) inner_down();
//color("Pink")  translate([300, 300, 2*z_cut-h1]) rotate([0, 180, 0]) inner_up();




// ============================================================================
// PIÈCE DE LIAISON VERTE — VERSION POUR INSERT CYLINDRIQUE DE SERRAGE
// ============================================================================
//du trou (4mm) et de sa profondeur (5mm).

// --- PARAMÈTRES DE L'INSERT MÉTALLIQUE ---
insert_diametre = 4;  // Diamètre extérieur de ton insert métallique
insert_hauteur  = 5;  // Longueur/Profondeur de l'insert dans le moyeu
insert_rayon     = insert_diametre / 2;
height_shaft = 24; // Longueur à ajuster pour englober le shaft 
diameter_shaft = 8;
diameter_link = diameter_shaft + 2*insert_hauteur;
echo("radius link: ", diameter_link/2);
fillet_r = 8;
// --- TOLERANCE SHAFT ---
tol_shaft = 0.5;

module fixation_piece_corrected() {
    // Le rayon externe s'aligne sur la collerette de montage des hoops
    flange_r = ushape_R + dy + ushape_oRings_dist/2; 
    difference() {
        // Disque de liaison principal
        cylinder(h=h1, r=flange_r, center=false, $fn=100);
        
        // Synchronisation avec les 6 trous M3 à la distance ushape_R
        for (i=[1:6]) {
            rotate([0, 0, i*60])
                translate([0, ushape_R, -1]) {
                    cylinder(h=h1 + 2, r=1.6, center=false, $fn=50);
                    translate([0, 0,  1])
                        mounting_hole_m3_tall();
                }
        }
    }
}

module conge() {
    translate([0, 0, -2*h1])
            difference() {
            rotate_extrude($fn=80)
            difference() {
                square([fillet_r+ diameter_link/2, fillet_r]);
                translate([fillet_r + diameter_link/2, 0]) circle(r=fillet_r, $fn=50);
            }
            cylinder(h=fillet_r, r=diameter_link/2);
            };
}

module final_link() {
    translate([0, 0, -h1])
    difference(){
        union() {
            // Collerette
            translate([0, 0, -height_shaft])
            cylinder(h=height_shaft+h1, r=diameter_link/2, center=false, $fn=80);
            // Disque 
            fixation_piece_corrected();
            // Congé
            conge();
            };

        // Trou du shaft
        translate([0, 0, -height_shaft])
        cylinder(h=height_shaft+h1, r=diameter_shaft/2+tol_shaft, center=false, $fn=50);
        
        // Trou de serrage
        translate([0, diameter_link/2, -holder_width/2])
        rotate([90, 0, 0])
        cylinder(h=insert_hauteur+3, r=insert_diametre/2, center=false, $fn=50);
    };    
}

// translate([0, (8+2*insert_hauteur)/2, -20])
//         rotate([90, 0, 0])
//         cylinder(h=insert_hauteur+3, r=insert_diametre/2, center=false, $fn=80);



echo("h1: ", h1);

//translate([0, 0, h1]) fixation_piece();
//translate([0, 0, -30+2*h1]) shaft();
//color("Red", 50) translate([0, 0, -30+2*h1]) center_link();

/// =============================================================================
// MOTOR SUPPORT WITH STIFFENERS
// ==================================================================

height_hole = 133.5;
holder_height = 180;
holder_length = 250;
holder_width = 20;
triangle_height = 95;
assembly_hole_radius = 5;

// FIXATIONS
width_fixation  = 20;
height_fixation = 60;
height_fixation_assembly_holes = 40;
length_fixation = 70;
// Position X centrée sur le bras pour tous les trous
center_x = (length_fixation - width_fixation) / 2;


// --- PARAMÈTRES DES NERVURES DE RIGIDITÉ ---
rib_thickness = 8;  // Épaisseur de la nervure
rib_height = 40;     // Hauteur de la nervure
rib_length = 120;    // Longueur de la nervure

module triangle() {
    // Augmentation légère de la hauteur pour un rendu propre (évite le Z-fighting)
    linear_extrude(height=holder_width + 2)
        polygon(points=[[0,0], [0,triangle_height],[triangle_height,triangle_height]]);
}

module support_box() {
    difference(){
        translate([-holder_length/2, -height_hole, -holder_width]) 
            cube([holder_length, holder_height, holder_width]);
        
        // Ajustement des Z pour que la soustraction soit parfaite
        translate([-holder_length/2, holder_height-(height_hole+triangle_height), -holder_width - 1]) 
            triangle();
        
        translate([holder_length/2-triangle_height,triangle_height/2, -holder_width - 1]) 
            rotate([0, 0, -90]) 
            triangle();
    }
}

module support_holes(){
    for (i=[1:4]){
        rotate([0,0,i*90])
            translate([23.57, 23.57, -1]) {
                cylinder(h=holder_width + 2, r=2.2, center=false, $fn=50);
                translate([0, 0, holder_width - 4 + 1])
                    cylinder(h=4.1, r=3.5, center=false, $fn=50);
        }   
    }
    // Trou central
    translate([0, 0, -1])
        cylinder(h=holder_width + 2, r=(diameter_link+2)/2, center=false,$fn=50);
}

// ============================================================
// TROU DE PASSAGE M6 (6.5 mm)
// ============================================================
module hole_M6() {
    cylinder(h=100, r=3.25, center=true, $fn=50);
}

module motor_support() {
    difference() {

        // La plaque de base avec ses perçages
        difference() {
            support_box();
            translate([0, 0, -holder_width]) support_holes();
        };

        scale(1.03) conge();

        translate([0, holder_height-height_hole, -holder_width/2 + 1])
        rotate([90, 0, 0])
        cylinder(h=holder_height-height_hole, r=5, center=false, $fn=50);

        translate([-100, -height_hole + height_fixation_assembly_holes, 0])
        hole_M6();

        translate([+100, -height_hole + height_fixation_assembly_holes, 0])
        hole_M6();
    }

}


// ============================================================
// BRAS DE FIXATION — équerre L reliant le support au sol + renfort
// ============================================================
module bras_fixation(z_trous, length_fixation = length_fixation) {

    difference() {
        translate([length_fixation/2, 0, 0])
        rotate([0, -90, 0])
        union() {
            // Face collée à la plaque
            cube([2*h1, height_fixation, width_fixation]);
            // Partie horizontale
            cube([length_fixation, 8, width_fixation]);
            // Renfort triangulaire
            difference() {
                cube([length_fixation/4, height_fixation/2, width_fixation]);
                translate([length_fixation/4+h1, height_fixation/3 + h1])
                    linear_extrude(height = width_fixation)
                        circle(r=17, $fn=80);
            }
        }
        translate([center_x, height_fixation_assembly_holes, 0])
        hole_M6();

        // Perçages M6 sur la grille 25 mm
        for (z = z_trous) {
            translate([center_x, h1/2, z])
                rotate([90, 0, 0])
                    hole_M6();
        }

    }
}


// ============================================================
// SUPPORT MOTEUR AVEC SES FIXATIONS AU SOL
// ============================================================

module fixations_with_holes() {
    L_fix          = length_fixation;          // allongé de 90 → 95 mm pour la marge au trou extérieur
    width_fixation = 20;

    center_x = (L_fix - width_fixation) / 2;   // = 37.5
    X_gauche = -100 - center_x;                 // trous X à -100 mm (!grille)
    X_droite =  100 - center_x;                 // trous X à +100 mm (!grille)

    union() {
        // --- Fixations arrière  (trous globaux : Z = -50 et -75 mm) ---
        rotate([0, 180, 0])
        translate([0, 0, holder_width]) {
            translate([X_gauche, -height_hole, -h1])
                bras_fixation([29, 54], L_fix);
            translate([X_droite, -height_hole, -h1])
                bras_fixation([29, 54], L_fix);
        }

        // --- Fixations avant   (trous globaux : Z = +50 et +75 mm) ---
        translate([X_gauche, -height_hole, -h1])
            bras_fixation([59-25, 84-25], L_fix);
        translate([X_droite, -height_hole, -h1])
            bras_fixation([59-25, 84-25], L_fix);
    }
}

module motor_support_with_fixations() {
    L_fix          = 70;          // allongé de 90 → 95 mm pour la marge au trou extérieur
    width_fixation = 20;

    center_x = (L_fix - width_fixation) / 2;   // = 37.5
    X_gauche = -100 - center_x;                 // trous X à -100 mm (!grille)
    X_droite =  100 - center_x;                 // trous X à +100 mm (!grille)

    translate([0, 0, -h1 - 1])
    union() {
        motor_support();

        fixations_with_holes();
    }
}


// ============================================================
// GRILLE DE RÉFÉRENCE — trous M6 pas 25 mm (visualisation)
// ============================================================
module grille() {
    translate([0, -2*h1-2, -height_hole])
    for (j = [-10:10]) {
        translate([0, j*25, 0])
        for (i = [-10:10]) {
            translate([i*25, 0, 0])
                hole_M6();
        }
    }
}

//grille();

// ==========================================
// MODULES OUTILS POUR L'ÉCLATÉ ET LA 2D
// ==========================================

module mise_en_plan_3_vues() {
    // Génère un triptyque de plans 2D orthogonaux style dessin industriel
    
    // 1. Vue de Face (Plan XY d'origine)
    translate([0, 0, 0]) 
        projection(cut = false) rotate([0, 0, 0]) children(0);
    
    // 2. Vue de Dessus (Plan XZ) - décalée vers le haut
    translate([0, 150, 0]) 
        projection(cut = false) rotate([90, 0, 0]) children(0);
    
    // 3. Vue de Profil (Plan YZ) - décalée vers la droite
    translate([200, 0, 0]) 
        projection(cut = false) rotate([0, 90, 0]) children(0);
}

module assembly(explode = 0) {
    dist_trous = 25;

    // 1. Le bloc moteur central bouge vers l'arrière 
    translate([0, 0, -explode]) 
        color("#F23E35") motor_support_with_fixations();

    // 2. Les anneaux (Outer Down et Up) s'écartent sur l'axe vertical Z
    translate([0, 0, explode]) {
        color("Cyan") outer_down();
        color("Cyan") outer_up();}

    // 3. Les anneaux (Inner Down et Up) s'écartent sur l'axe vertical Z
    translate([0, 0, 2*explode]){
        color("Orange") inner_down();
        color("Orange") inner_up();}

    //4. Le link ne bouge pas
    color("Green") final_link();
}

module orienter_axes() {
    rotate([90, 0, 0]) children();
}
