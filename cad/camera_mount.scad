// ============================================================================
//  SUPPORT CAMERA A RAIL — Ball-in-Double-Hoop
// ----------------------------------------------------------------------------
//
//  CONCEPTION
//    [1] Liaison prismatique a 1 DDL (queue d'aronde) : le profil trapezoidal
//        contraint SIMULTANEMENT le soulevement et le mouvement lateral.
//    [2] Angle : 60 deg (flancs a 30 deg de la
//        verticale) : valeur normalisee pour les PETITES glissieres.
//    [3] Alignement garanti PAR CONSTRUCTION : 4 boulons M6 symetriques sur les
//        colonnes X = +/-25 mm de la grille  ->  centrage sur X = 0 et
//        parallelisme a l'axe du rotor.
//    [4] Hauteur d'axe optique = 133.5 mm = niveau de l'axe du rotor
//        (height_hole dans SCAD de l'assemblage principal).
//
//  A VERIFIER avant impression (marqué "!VERIF") :
//    - position exacte du connecteur CSI sur la carte (marge estimee, cf. cable_channel)
//    - jeu fonctionnel de la queue d'aronde (dt_cl)?
//
//  SOURCE DIMENSIONS CAMERA (mise a jour) :
//    Dessin mecanique officiel "RASPBERRY PI CAMERA MODULE V2.1",
//    ref. RPI-CAM-V2_1, dessine par Mike Stimson, approuve James Adams,
//    Raspberry Pi (Trading) Ltd., 12/11/2015.
//    - Carte : 25 x 23.862 mm, coins arrondis R2mm, epaisseur ~1.5mm
//    - 4x trous de fixation, diametre 2.2mm
//    - Entraxe trous : 21mm (axe X, cam_mount_dx) x 12.5mm (axe Z, cam_mount_dz)
//    - Le motif de trous N'EST PAS centre sur l'axe optique dans la direction
//      Z (23.862mm) : distance trou-cote-cable -> centre optique = 3.038mm ;
//      distance trou-cote-oppose -> centre optique = 9.462mm (3.038+9.462=12.5,
//      coherent avec l'entraxe). D'ou l'offset cam_hole_off_z ci-dessous.
// ============================================================================

$fn = 48;

// ---------------------------------------------------------------------------
//  PARAMETRES GENERAUX
// ---------------------------------------------------------------------------

// --- Grille aluminium (source : SCAD assemblage principal) ---
grid_pitch   = 25;     // [mm] pas de la grille de trous M6
m6_clear     = 6.5;    // [mm] diam trou de passage M6 (jeu)   (= hole_M6, r=3.25)
m6_head_dia  = 11;     // [mm] diam lamage tete cylindrique M6
m6_head_h    = 6;      // [mm] hauteur tete M6 (CHC)

// --- Alignement optique ---
axis_height  = 133.5;  // [mm] hauteur du centre optique / face grille   [4]

// --- RaspiCam V2 --- (dessin mecanique officiel RPI-CAM-V2_1, Raspberry Pi, 2015)
cam_w          = 25;      // [mm] largeur carte (X)
cam_h          = 23.862;  // [mm] hauteur carte (Z)                          // <-- MODIFIE (etait 24, valeur exacte du dessin)
cam_t          = 1.5;     // [mm] epaisseur carte
cam_cl         = 0.4;     // [mm] jeu camera
cam_mount_dx   = 21;      // [mm] entraxe trous fixation (horizontal)   — CONFIRME dessin officiel
cam_mount_dz   = 12.5;    // [mm] entraxe trous fixation (vertical)     — CONFIRME dessin officiel
cam_mount_d    = 2.2;     // [mm] diam trou fixation                    — CONFIRME dessin officiel (etait suppose M2 autotaraudeur ~2.2mm, coherent)
cam_hole_off_z = 3.212;   // [mm] decalage du CENTRE du motif de trous par rapport au centre   // <-- NOUVEAU
                           //      optique, vers le cote OPPOSE au cable CSI (dessin off. :
                           //      9.462mm trou-haut<->objectif, 12.5-9.462=3.038mm trou-bas(cable)<->objectif
                           //      => offset = (9.462-3.038)/2 = 3.212mm)
lens_off_z     = 0;       // [mm] decalage objectif / centre carte (X/Y, non concerne par cam_hole_off_z)

// --- DoveTail ---
dt_angle   = 60;                          // [deg] angle inclus normalise   [2]
dt_flank   = 90 - dt_angle;               // [deg] flanc / verticale = 30 deg
dt_h       = 7;                           // [mm] hauteur de la queue d'aronde
dt_narrow  = 16;                          // [mm] largeur cote etroit (base)
dt_wide    = dt_narrow + 2*dt_h*tan(dt_flank);  // [mm] largeur cote large (haut)
dt_cl      = 0.35;                        // [mm] jeu fonctionnel femelle (FDM) — valeur nominale
dt_groove_w = 4;                          // [mm] largeur rainure superieure (vis de blocage)
dt_groove_d = 1.5;                        // [mm] profondeur rainure superieure

// --- Rail ---
rail_base_w = 64;      // [mm] largeur base (X) : accueille boulons +/-25 + marge
rail_base_t = 8;       // [mm] epaisseur base (Z) : permet le lamage tete M6
rail_len    = 125;     // [mm] longueur (Y) le long de l'axe optique
bolt_dx     = 25;      // [mm] boulons sur colonnes X = +/-25 (grille)     [3]
bolt_dy     = 50;      // [mm] boulons sur lignes  Y = +/-75 (span 150 = 6 pas)

// --- Chariot ---
blk_w      = 36;       // [mm] largeur bloc queue d'aronde (X)
blk_h      = 14;       // [mm] hauteur bloc (Z) au-dessus de la base rail
blk_y0     = -30;      // [mm] bloc : bord avant (Y)
blk_y1     =  18;      // [mm] bloc : bord arriere (Y)   -> longueur 45 mm
head_x     = 32;       // [mm] largeur tete/colonne (X)
head_y0    = 0;        // [mm] face AVANT de la tete (plan d'appui carte)
head_y1    = 18;       // [mm] face ARRIERE de la tete
boss_x     = 24;       // [mm] extension laterale du bossage vis de blocage
set_insert_d = 4;      // [mm] diam insert laiton (vis M3 de blocage / gib)   [1]
set_insert_h = 5;      // [mm] profondeur insert

// --- Cablage (ruban CSI cache) ---
cable_w    = 18;       // [mm] largeur du canal (ruban CSI ~16 mm)
cable_d    = 8;        // [mm] profondeur du canal (Y)

// --- Zone haute (tete) ---
head_top   = axis_height + cam_h/2 + 4;   // [mm] sommet de la tete
blk_top    = rail_base_t + blk_h;         // [mm] Z du dessus du bloc dovetail

// ---------------------------------------------------------------------------
//  ECHOS DE CONTROLE
// ---------------------------------------------------------------------------
car_len   = blk_y1 - blk_y0;                 // longueur chariot
travel    = rail_len - car_len;              // course utile (Y)
echo(str("Largeur queue d'aronde (haut) dt_wide = ", dt_wide, " mm"));
echo(str("Longueur chariot = ", car_len, " mm ; course utile = ", travel, " mm"));
echo(str("Placer le rail sur la grille pour que car_y=0 donne d ~ 305 mm (milieu)"));
echo(str("Hauteur camera = ", axis_height, "mm"));
echo(str("Centre motif de trous = ", axis_height + cam_hole_off_z, " mm (offset ", cam_hole_off_z, " mm / axe optique)"));  // <-- NOUVEAU

// ===========================================================================
//  MODULES DE BASE
// ===========================================================================

// Profil 2D de la queue d'aronde (male), base en y=0, largeur suivant x.
// off > 0 : offset uniforme (jeu) pour generer la rainure femelle.

module dovetail_2d_groove(off = 0) {
    nn = dt_narrow/2;
    ww = dt_wide/2;
    gw = dt_groove_w/2;
    offset(delta = off)
        polygon([
            [-nn, 0], [nn, 0],
            [ww, dt_h],
            [gw, dt_h], [gw, dt_h - dt_groove_d],
            [-gw, dt_h - dt_groove_d], [-gw, dt_h],
            [-ww, dt_h]
        ]);
}

module dovetail_2d(off = 0) {
    nn = dt_narrow/2;
    ww = dt_wide/2;
    gw = dt_groove_w/2;
    offset(delta = off)
        polygon([
            [-nn, 0], [nn, 0],
            [ww, dt_h],
            [gw, dt_h], [-gw, dt_h],
            [-ww, dt_h]
        ]);
}

// Solide queue d'aronde : longueur le long de Y, base a z=0.
//   rotate([90,0,0]) envoie l'axe d'extrusion (Z) sur -Y et la hauteur (y) sur Z.
module dovetail_solid_groove(length, off = 0) {
    rotate([90, 0, 0])
        linear_extrude(height = length, center = true)
            dovetail_2d_groove(off);
}

module dovetail_solid(length, off = 0) {
    rotate([90, 0, 0])
        linear_extrude(height = length, center = true)
            dovetail_2d(off);
}

// Trou de passage M6 + lamage tete (perce selon Z, tete noyee par le haut).
module m6_hole() {
    translate([0, 0, -1])
        cylinder(h = rail_base_t + 2, d = m6_clear);
    translate([0, 0, rail_base_t - m6_head_h])
        cylinder(h = m6_head_h + 1, d = m6_head_dia);
}

// ===========================================================================
//  RAIL
// ===========================================================================
module rail() {
    difference() {
        union() {
            // Base plate posee sur la grille
            translate([-rail_base_w/2, -rail_len/2, 0])
                cube([rail_base_w, rail_len, rail_base_t]);
            // Queue d'aronde sur toute la longueur
            translate([0, 0, rail_base_t])
                dovetail_solid_groove(rail_len);
        }
        // 4 boulons M6 aux coins, sur la grille (X = +/-25, Y = +/-75)  [3]
        for (sx = [-1, 1], sy = [-1, 1])
            translate([sx*bolt_dx, sy*bolt_dy, 0]) m6_hole();
    }
}

// ===========================================================================
//  CHARIOT PORTE-CAMERA
// ===========================================================================

// -- 1. Bloc a queue d'aronde (partie qui coulisse) --
module dovetail_block() {
    translate([-blk_w/2, blk_y0, rail_base_t])
        cube([blk_w, car_len, blk_h]);
}

// Rainure femelle = dovetail male + jeu, traversante en Y (coulissement).
// Ouverte vers le bas (offset descend sous la face) pour un montage franc.
// off    = jeu fonctionnel (parametrable pour les eprouvettes de calibration,
//          defaut = dt_cl pour ne pas casser les appels existants de carriage()).
// length = longueur de la coupe (defaut = car_len + 20, comme a l'origine).
module dovetail_slot(off = dt_cl, length = car_len + 20) {
    translate([0, 0, rail_base_t - 0.01])
        dovetail_solid(length, off);
}

// -- 2. Tete (colonne creuse + plaque porte-camera) --
module head_solid() {
    translate([-head_x/2, head_y0, blk_top])
        cube([head_x, head_y1 - head_y0, head_top - blk_top]);
}

// -- 3. Nervure arriere (rigidite fore/aft ; charge camera ~3 g negligeable,
//        objectif = raideur anti-vibration, pas resistance) --
module gusset() {
    gh = 42;  gl = 12;  th = 8;
    translate([-th/2, head_y1, blk_top])
        rotate([0, 0, 90]) rotate([90, 0, 0])
            linear_extrude(height = th)
                polygon([[0, 0], [gl, 0], [0, gh]]);
}

// -- 4. Bossage lateral pour la vis de blocage (gib) --
module set_screw_boss() {
    translate([blk_w/2 - 0.01, 0, rail_base_t + 1])
        cube([boss_x - blk_w/2, 12, dt_h + 2]);
}

// -- Soustractions de la tete --

// Logement / appui de la carte : leger evidement de reperage (0.8 mm)
// Centre sur l'axe optique (axis_height + lens_off_z) — INCHANGE, c'est l'objectif
// qui doit s'aligner sur l'axe du rotor, independamment de la position du motif de trous.
module camera_pocket() {
    translate([-(cam_w + cam_cl)/2, -0.01, axis_height + lens_off_z - (cam_h + cam_cl)/2])
        cube([cam_w + cam_cl, 1 + 0.01, cam_h + cam_cl]);
}

// 4 trous de fixation M2 de la carte, perces dans la face avant (Y=0).
// Motif recentre sur (axis_height + cam_hole_off_z) : le motif de trous n'est
// PAS centre sur l'axe optique (cf. dessin mecanique officiel, en tete de fichier).
module camera_mount_holes() {                                                    // <-- MODIFIE
    for (sx = [-1, 1], sz = [-1, 1])
        translate([sx*cam_mount_dx/2, -1, axis_height + cam_hole_off_z + sz*cam_mount_dz/2])  // <-- MODIFIE (+ cam_hole_off_z)
            rotate([-90, 0, 0]) cylinder(h = 9, d = cam_mount_d);
}

// Chemin de cable CSI cache : evidement connecteur (avant, bas de carte)
// -> canal vertical interne -> sortie arriere en bas (vers le Raspberry Pi).
// Le connecteur est du cote du trou le plus proche de l'axe optique (3.038mm,
// cf. dessin officiel) ; z_conn recale en consequence (marge de 4mm vers le
// bord de carte encore estimee, non cotee sur le dessin -> !VERIF).
module cable_channel() {
    z_conn = axis_height - 3.038 - 4;   // proche du trou cote cable ; marge 4mm estimee  // <-- MODIFIE  !VERIF
    // (a) evidement connecteur, ouvert sur la face avant
    translate([-cable_w/2, -0.01, z_conn - 3])
        cube([cable_w, 5, 10]);
    // (b) canal vertical interne (relie l'evidement au bas de la tete)
    translate([-cable_w/2, 4, blk_top - 0.01])
        cube([cable_w, cable_d, z_conn - blk_top + 6]);
    // (c) sortie arriere en bas
    translate([-cable_w/2, head_y1 - 6, blk_top - 0.01])
        cube([cable_w, 7, 12]);
    // (d) evidement connecteur sous les vis du bas
    translate([-cam_w/2, -0.01, z_conn - 8])
        cube([cam_w, 6, 10]);
}

// Vis de blocage M3 : percage VERTICAL depuis le dessus du chariot
module set_screw_hole() {
    z_bottom = rail_base_t + dt_h - dt_groove_d;
    translate([0, 5, z_bottom])
        cylinder(h = blk_top - z_bottom + 2, d = set_insert_d);
}

module carriage(car_y = 0) {
    translate([0, car_y, 0])
    difference() {
        union() {
            dovetail_block();
            head_solid();
        }
        dovetail_slot();
        camera_pocket();
        camera_mount_holes();
        cable_channel();
        translate([0, -15, 0]) 
        set_screw_hole();
    }
}

// ===========================================================================
//  EPROUVETTES DE CALIBRATION DU JEU FONCTIONNEL (queue d'aronde)
// ---------------------------------------------------------------------------
//  But : imprimer plusieurs blocs femelles (= dovetail_block() + rainure)
//  a des jeux differents, pour les faire glisser sur le rail deja imprime
//  et determiner le jeu FDM fonctionnel avant d'imprimer la tete complete.
//  Valeurs testees : encadrent la valeur nominale actuelle (dt_cl = 0.35 mm).
// ===========================================================================

dt_cl_test_values = [0.20, 0.35, 0.50];   // [mm] jeu serre / nominal / lache
test_block_len     = 20;                  // [mm] longueur des eprouvettes (au lieu de car_len = 45 mm)

// Bloc court dedie aux eprouvettes — independant de dovetail_block()/car_len,
// centre en Y=0 (peu importe sa position exacte le long du rail pour ce test).
module test_block(length = test_block_len) {
    translate([-blk_w/2, -length/2, rail_base_t])
        cube([blk_w, length, blk_h]);
}

// Module generique parametre par le jeu — utilise par base_jeu1/2/3.
module base_jeu(off, length = test_block_len) {
    difference() {
        test_block(length);
        dovetail_slot(off, length + 10);
    }
}

module base_jeu1() { base_jeu(dt_cl_test_values[0]); }  // jeu serre
module base_jeu2() { base_jeu(dt_cl_test_values[1]); }  // jeu nominal (= dt_cl)
module base_jeu3() { base_jeu(dt_cl_test_values[2]); }  // jeu lache

// Etiquette gravee (en creux, face superieure) pour identifier chaque
// eprouvette apres impression — sinon impossible de les distinguer une fois
// sorties de l'imprimante.
module base_jeu_label(off, length = test_block_len) {
    translate([0, length/2 - 6, blk_top - 0.4])
        linear_extrude(height = 0.6)
            text(str(off), size = 5, halign = "center", valign = "center");
}

module base_jeu_labeled(off, length = test_block_len) {
    difference() {
        base_jeu(off, length);
        base_jeu_label(off, length);
    }
}

// Les 3 eprouvettes posees cote a cote sur le plateau, pretes pour un seul job
// d'impression (espacement = largeur du bloc + marge).
module base_jeu_all() {
    spacing = blk_w + 10;
    for (i = [0:2])
        translate([(i - 1) * spacing, 0, 0])
            base_jeu_labeled(dt_cl_test_values[i]);
}

// ===========================================================================
//  GRILLE DE REFERENCE (visualisation de l'alignement)
// ===========================================================================
module grid_ref() {
    color([0.6, 0.6, 0.6, 0.25])
    translate([0, 0, -2])
    for (i = [-2 : 2], j = [-4 : 4])
        translate([i*grid_pitch, j*grid_pitch, 0])
            cylinder(h = 2, d = m6_clear);
}

