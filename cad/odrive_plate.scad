// ============================================================================
//  PLAQUE ODrive S1 — version simple
// ============================================================================
//  Une plaque, taille du driver + marge pour 4 trous M6 (fixation plaque alu).
//  Sur la plaque : 4 petits trous pour visser le driver directement.
//
//  Cotes des trous M3 confirmées sur le modèle Onshape officiel ODrive S1
//  (mesure directe, voir commentaire dans la section "Trous M3" ci-dessous).
// ============================================================================

/* ---- Paramètres ---- */

// Empreinte PCB [mm] — source: page produit ODrive S1 (51x64mm)
pcb_L = 64.0;
pcb_W = 51.0;

// Marge ajoutée autour du PCB pour loger les trous M6 (calculée plus bas
// à partir de l'espacement M6 souhaité, ne pas éditer directement)

// Épaisseur de la plaque
thick = 4.0;

// --- Trous M6 vers la plaque aluminium ---
// La grille de la plaque alu a un pas de 25mm en X et en Y (indépendamment).
// Pour tomber pile sur des trous existants quelle que soit la position de
// la plaque sur la grille, l'espacement centre-à-centre des 4 trous M6 doit
// être un multiple de 25mm SUR CHAQUE AXE (X et Y peuvent avoir des
// multiples différents, la grille n'a pas besoin d'être carrée).
m6_hole_d    = 6.5;
m6_spacing_x = 75.0;   // 3 x 25mm
m6_spacing_y = 50.0;   // 2 x 25mm
m6_edge_clear = 8.0;   // marge mini entre trou M6 et bord de plaque [mm]

// --- Trous M3 de fixation du driver (sur le PCB) ---
// Cotes mesurées sur le modèle Onshape officiel (outil Measure) : pattern rectangulaire 60.0mm (X) x 44.0mm (Y), centré sur l'empreinte 64x51mm
// -> inset 2.0mm (X) / 3.5mm (Y) depuis chaque bord. Diamètre mesuré: 3.4mm.
pcb_hole_d = 3.4;   // jeu M3, mesuré (Diameter: 3.400mm sur Onshape)
pcb_holes = [
    [2.0,        3.5],          // coin 1
    [pcb_L-2.0,  3.5],          // coin 2
    [2.0,        pcb_W-3.5],    // coin 3
    [pcb_L-2.0,  pcb_W-3.5],    // coin 4
];

/* ---- Géométrie ---- */

// La plaque doit contenir le pattern M6 (spacing + marge de chaque côté)
// ET le PCB (avec un minimum de marge autour) — on prend le plus grand des deux.
plate_L = max(m6_spacing_x + 2*m6_edge_clear, pcb_L + 10);
plate_W = max(m6_spacing_y + 2*m6_edge_clear, pcb_W + 10);

// PCB centré dans la plaque
off_x = (plate_L - pcb_L) / 2;
off_y = (plate_W - pcb_W) / 2;

// Trous M6 centrés dans la plaque, espacés de m6_spacing_x/y
m6_center_x = plate_L / 2;
m6_center_y = plate_W / 2;

module plate() {
    union() {
        difference() {
            cube([plate_L, plate_W, thick]);

            // Trous M6, espacés de m6_spacing_x (X) et m6_spacing_y (Y), centrés
            m6_positions = [
                [m6_center_x - m6_spacing_x/2,  m6_center_y - m6_spacing_y/2],
                [m6_center_x + m6_spacing_x/2,  m6_center_y - m6_spacing_y/2],
                [m6_center_x - m6_spacing_x/2,  m6_center_y + m6_spacing_y/2],
                [m6_center_x + m6_spacing_x/2,  m6_center_y + m6_spacing_y/2],
            ];
            for (p = m6_positions)
                translate([p[0], p[1], -1])
                    cylinder(d = m6_hole_d, h = thick + 2, $fn = 32);
            }
        // Trous M3 (positions du PCB, décalées de l'offset)
        translate([0, 0 , thick])
        for (h = pcb_holes)
            translate([off_x + h[0], off_y + h[1], -1])
                union() {
                    difference() {
                        union() {
                        cylinder(d = pcb_hole_d, h = thick+2, $fn = 32);
                         translate([0, 0, thick+pcb_hole_d/2]) sphere(d=pcb_hole_d+0.5, $fn=50);
                        }
                        translate([-pcb_hole_d/2-2, -pcb_hole_d/8,0])
                        cube([pcb_hole_d+4, pcb_hole_d/4, thick+8]);
                    }
                   
                }
            
    }
}

plate();
