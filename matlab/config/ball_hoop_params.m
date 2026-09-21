function params = ball_hoop_params()
%BALL_HOOP_PARAMS  Physical parameters of the ball-in-double-hoop system.
%
%   params = ball_hoop_params()
%
%   Inputs:  none
%   Outputs: params (struct) — all physical quantities in SI units
%
%   Contains ONLY physics: geometry, masses, inertias, friction coefficients.
%   No numerical settings (those belong in sim_config), no controller gains
%   (those belong in controller structs), no input signals (those belong in
%   scenario files).
%
%   To swap the physical system replace this file only. Every downstream
%   module reads exclusively from the params struct — no hard-coded values
%   anywhere else.
%
%   PROVENANCE. Every value below is tagged with how it was obtained, using
%   the author's own parameter table for the as-built prototype:
%     cst      — physical constant or manufacturer spec
%     measure  — directly measured on the prototype
%     model    — from the CAD/design geometry
%     formula  — derived from other entries in this file
%     exp.     — identified experimentally (bench test)
%   ball_friction was the last PROVISIONAL entry here; experiment E2 has
%   since identified it, so every value below now has a stated provenance.
%
%   Reference geometry: Zemanek (2017) ball-in-double-hoop, re-measured for
%   this prototype.

    %% --- Environment ---
    params.gravity = 9.81;                                  % [m/s^2]   cst

    %% --- Ball ---
    params.ball_mass    = 67e-3;                            % [kg]  measure
    params.ball_radius  = 12.5e-3;                          % [m]   measure

    % Solid-sphere formula; valid because the ball has uniform density.
    params.ball_inertia = (2/5) * params.ball_mass * params.ball_radius^2;
                                                            % [kg.m^2]  formula

    % Identified in experiment E2 (free-oscillation decay of the ball in the
    % outer hoop), 6.484e-6 +/- 0.84e-6, superseding the earlier 2.57e-6
    % placeholder. It sets the decay envelope, not the frequency, so it is
    % the parameter a free-oscillation fit constrains best.
    %
    % E2 fitted this against the AS-BUILT two-rail geometry, and that is
    % the geometry hoop_geometry.m now applies (r_roll = 10.639 mm,
    % R_eff = 91.58 mm on the outer track), so b_b and the geometry it is
    % paired with agree. The linearised open loop at psi = 0 comes out at
    % zeta = 0.0332, E2's own decrement.
    params.ball_friction = 6.484e-6;                        % [N.m.s/rad]  exp. (E2)

    %% --- O-ring tracks (two-rail contact geometry) ---
    % The ball is not carried by a cylindrical surface: each track is a pair
    % of O-rings, axially separated, so the ball sits in a vee and touches
    % each track at TWO points symmetric about its mid-plane. Two
    % consequences, both stiffening the model (thesis eq. 78-79):
    %
    %   Delta  = sqrt((Rb + rc)^2 - (h/2)^2)      radial ball-centre offset
    %   R_eff  = Rg -/+ Delta                     (- concave, + convex)
    %   r_roll = Rb * Delta / (Rb + rc)           rolling radius, < Rb
    %
    % with Rg the O-ring CENTRELINE radius of the track. For the as-built
    % prototype this gives Delta = 11.92 mm, r_roll = 10.64 mm, and
    % R_eff = 91.6 / 67.4 / 32.6 mm for the outer, inner-outside and
    % inner-inside tracks. The ball rolls on a chord of itself rather than a
    % diameter, so it must spin faster to travel the same distance: its
    % effective rotational inertia rises by the factor
    % k = 1 + (2/5)(Rb/r_roll)^2 = 1.552, against 1.4 for the single-contact
    % idealisation.
    %
    % STATUS: applied. hoop_geometry.m derives Delta and r_roll from the two
    % measurements below and returns the two-rail R_eff and r_roll;
    % rolling_matrices.m divides every no-slip ratio by that r_roll, keeping
    % ball_radius only for the ball's own inertia (2/5)m*Rb^2 and for the
    % gravity term. There is deliberately NO params.ball_rolling_radius
    % field: r_roll is a consequence of rc and h, and storing it separately
    % would be a second source of truth free to drift out of step with them.
    %
    % Do not take the r_roll consequence without the Delta one. The two come
    % from the same tangency construction, and they push the natural
    % frequency in opposite directions -- k rises 11 %, R_eff falls 8.4 %,
    % and they very nearly cancel: 1.3323 Hz single-contact against
    % 1.3221 Hz two-rail, a 0.8 % shift, not the ~5 % that changing r_roll
    % alone would suggest. The bench measured 1.3144 Hz (E2).
    params.oring_cord_radius   = 1.5e-3;                    % [m]   measure, cord section rc
    params.oring_axial_spacing = 14.7e-3;                   % [m]   measure, separation h
                                                            %       (NOT axial_spacing below,
                                                            %        which separates hoop planes)

    % Static friction coefficient of the steel-ball / O-ring contact. It
    % sets the no-slip bound on hoop acceleration (contact_forces.m,
    % analysis/verification/slip_acceleration_limit.m), which is the binding
    % actuator constraint for T3/T4 once the second supply raises tau_max
    % above ~0.15 N.m -- below that the torque limit binds first and this
    % value does not matter.
    %
    % NOT measured on this prototype. It is the value the bench code
    % already assumes (Codes/Python/t1_t2/t1_config_v3.py, U_SLIP_MAX),
    % carried here so the two sides agree rather than each hard-coding a
    % number. An inclined-plane test -- ball on a length of the same
    % O-ring stock, tilt until it slides, mu = tan(angle) -- settles it
    % directly.
    %
    % LITERATURE CHECK (2026-08-10), since the measurement is not
    % available: published dry elastomer-on-steel friction is HIGHER than
    % this, not lower. RoyMech's compilation gives static "Solids on
    % Rubber" as 1.0-4.0; general NBR figures quoted by suppliers sit in
    % 0.4-1.0; a reported test of uncoated NBR rings against a steel
    % counterface gives 0.48 at the lowest sliding speed rising to 1.51 at
    % 99 mm/s. Rubber friction is adhesion-dominated and routinely exceeds
    % 1, unlike the metal-on-metal intuition.
    %
    % The value is therefore KEPT AT 0.5 and re-labelled: it is a
    % conservative lower bound, not a best estimate. Raising it toward the
    % literature centre (~1.0) would relax the no-slip constraint and make
    % every T3/T4 result look better, which is the wrong direction to move
    % an unmeasured parameter. Results that depend on it are reported with
    % a sensitivity sweep over mu in [0.3, 0.6], never at this value alone
    % -- and the T3 looping is feasible across that whole range
    % (analysis/studies/T3/t3c_offline_study.m: 5.24 A at mu = 0.3 against 4.08 A
    % at 0.5), so the conclusions do not rest on it.
    params.ball_track_friction_coeff = 0.5;                 % [-]   CONSERVATIVE, not measured

    %% --- Outer hoop ---
    params.outer_hoop_radius         = 102.5e-3;            % [m]   model, inner (rolling) surface
    params.outer_hoop_inner_radius   = params.outer_hoop_radius;

    % O-ring centreline radius of the outer track: the rail sits proud of
    % the hoop's inner surface, so this is NOT outer_hoop_radius - rc and
    % cannot be derived from the values above -- it is measured as-built.
    % Its inner-hoop counterparts need no separate field: the inner hoop's
    % rails sit on its rod surfaces, so inner_hoop_outer_radius (55.5 mm)
    % and inner_hoop_inner_radius (44.5 mm) already ARE the centreline
    % radii of the convex and concave inner tracks.
    params.outer_hoop_oring_radius   = 103.5e-3;            % [m]   measure, Rg

    % Display only: outer_hoop_outer_radius sets the drawn hoop edge in
    % plot_trajectory.m / animate_system.m and enters no equation of motion.
    % Not in the prototype parameter table -- kept at its previous value.
    params.outer_hoop_thickness      = 10e-3;               % [m]
    params.outer_hoop_outer_radius   = params.outer_hoop_radius + params.outer_hoop_thickness;

    %% --- Inner hoop ---
    % inner_hoop_thickness is the FULL cross-section thickness (rod
    % diameter); the rolling surfaces sit at +/- half of it from the
    % mid-plane radius. This reproduces the prototype table's
    % "Ri,in/out = Ri -/+ t_i" with its t_i = 5.5 mm half-thickness:
    % 50 -/+ 5.5 = 44.5 / 55.5 mm.
    params.inner_hoop_radius         = 50e-3;               % [m]   model, mid-plane radius
    params.inner_hoop_thickness      = 11e-3;               % [m]   model
    params.inner_hoop_outer_radius   = params.inner_hoop_radius + params.inner_hoop_thickness/2;
    params.inner_hoop_inner_radius   = params.inner_hoop_radius - params.inner_hoop_thickness/2;

    %% --- Inner-hoop gap ---
    % The gap lets the ball transition between rolling_in_outside and
    % rolling_in_inside, and is the only path to free_fall from the inner hoop.
    params.hole_angular_width = pi/2;                       % [rad]  model, 90-degree opening
    params.hole_center_angle  = pi;                         % [rad]  centered at the bottom

    %% --- Hoop-pair geometry ---
    % Neither quantity is read anywhere else in the codebase, and neither
    % appears in the prototype parameter table; kept at their previous
    % values pending confirmation against the as-built assembly.
    params.axial_spacing = 30e-3;                           % [m]   gap between hoop planes
    params.hoop_width    = 30e-3;                           % [m]   hoop cross-section width

    %% --- Hoop assembly + motor inertia ---
    % Both inertias are bench-identified for the assembled prototype, which
    % supersedes the previous thin-ring estimate built from assumed hoop
    % masses (that formula, and the outer_hoop_mass / inner_hoop_mass values
    % it needed, are gone: a measured total beats a modelled one, and the
    % masses had no other reader).
    params.hoop_inertia        = 1.08e-3;                   % [kg.m^2]  exp.  (hoops alone)
    params.motor_rotor_inertia = 3.93e-5;                   % [kg.m^2]  exp.

    % Rotor inertia is lumped in with the hoops because the motor is rigidly
    % coupled to them; reproduces the prototype table's J_tot = 1.12e-3.
    params.hoop_motor_inertia = params.hoop_inertia + params.motor_rotor_inertia;
                                                            % [kg.m^2]  formula

    %% --- Motor ---
    params.motor_friction = 6.19e-5;                        % [N.m.s/rad]  exp., viscous

    % Coulomb (dry) friction torque, bench-identified. NOT currently wired
    % into the equations of motion: rolling_matrices.m and
    % free_fall_matrices.m model viscous damping only, and docs/VALIDATION.md
    % sec. 1.1 deliberately keeps Coulomb friction in the per-scenario
    % disturbance budget rather than the dynamics. Adding it is a model
    % change, not a parameter import -- it introduces a sign() discontinuity
    % at zero velocity that the stiff ODE solver and the CasADi planner
    % (which needs smooth derivatives) would both have to be adapted for.
    % Recorded here so the measured value is not lost.
    params.motor_coulomb_friction = 8.90e-3;                % [N.m]  exp.

    % Electrical characteristics. Not read by the mechanical simulation --
    % torque enters through cfg.tau_max (sim_config.m), which is itself
    % derived from torque_constant x the configured ODrive current limit.
    % Recorded here as the single source of truth for that derivation.
    params.motor_pole_pairs      = 7;                       % [-]        cst
    params.motor_kv              = 270;                     % [rpm/V]    cst
    params.motor_torque_constant = 0.031;                   % [N.m/A]    cst
    params.motor_phase_resistance = 38.487e-3;              % [Ohm]  exp. (datasheet: 39e-3)
    params.motor_phase_inductance = 16.40275e-6;            % [H]    exp. (datasheet: 16e-6)

    % Inertia the emulated ODrive inner velocity loop uses to convert the
    % acceleration setpoint u into a torque (see motor_inner_loop.m). The
    % ideal inner-loop model deliberately ignores the ball's coupling
    % through the rolling contact -- "just as the ODrive's velocity loop
    % would (it does not know about the ball)" -- so this is the hoop+rotor
    % inertia alone, not the mode-dependent M(1,1) rolling_matrices.m
    % computes while the ball is in contact. See docs/AUDIT.md item 1.5.
    params.total_inertia = params.hoop_motor_inertia;       % [kg.m^2]
end
