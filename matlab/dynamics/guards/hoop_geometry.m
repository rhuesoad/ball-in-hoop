function [surface, R_eff, R_track, coupling_sign, psi_eq, r_roll] = hoop_geometry(mode, params)
%HOOP_GEOMETRY  Single source of truth for the mode -> surface -> geometry
%               mapping used throughout the ball-hoop codebase.
%
%   [surface, R_eff, R_track, coupling_sign] = hoop_geometry(mode, params)
%   [surface, R_eff, R_track, coupling_sign, psi_eq, r_roll] = hoop_geometry(mode, params)
%
%   TWO-RAIL CONTACT. Each track is a pair of O-rings axially separated by
%   h, so the ball sits in a vee and touches the track at two points rather
%   than one (thesis eq. 78-79, and the derivation in ball_hoop_params.m's
%   O-ring section). Two quantities follow, and they are NOT the naive
%   R_track -/+ ball_radius and ball_radius that this function used to
%   return:
%
%       Delta  = sqrt((Rb + rc)^2 - (h/2)^2)
%       R_eff  = R_track -/+ Delta        (- concave, + convex)
%       r_roll = Rb * Delta / (Rb + rc)   < Rb
%
%   The two contact points define a rolling axis that is a chord of the
%   ball, not a diameter, so the ball must spin faster to cover the same
%   distance and its effective rotational inertia rises. Every no-slip
%   ratio downstream therefore divides by r_roll, never by ball_radius --
%   only the ball's own inertia (2/5)m*Rb^2 still uses Rb.
%
%   Inputs
%   ------
%   mode   : dynamic/control mode string --
%              'rolling_out' | 'rolling_in_outside' | 'rolling_in_inside'
%   params : physical parameters struct (see ball_hoop_params)
%
%   Outputs
%   -------
%   surface       : surface name used internally by rolling_matrices.m
%                   ('outer' | 'inner_outside' | 'inner_inside')
%   R_eff         : ball-centre orbit radius                          [m]
%   R_track       : O-ring CENTRELINE radius of this track [m]. Formerly
%                   named R_contact and equal to the hoop's rolling
%                   surface radius; with two rails the ball never touches
%                   that surface, so the quantity that sets the no-slip
%                   ratio is the rail centreline instead.
%   coupling_sign : sign of the M(1,2)/M(2,1) inertial-coupling term
%                   (-1 concave contact, +1 convex contact)
%   psi_eq        : (optional, 5th output) natural equilibrium angle
%                   for this mode [rad]. Not part of "hoop geometry"
%                   proper (it's a control-design quantity, the
%                   linearisation point each mode is naturally stable
%                   or unstable around) but kept here rather than
%                   reintroducing a second mode-keyed switch elsewhere
%                   for the sake of one field; the only caller that
%                   needs it (lqr_design.m) requests it explicitly.
%   r_roll        : (optional, 6th output) rolling radius [m] -- the
%                   perpendicular distance from a contact point to the
%                   ball's spin axis. Same for all three tracks, since it
%                   depends only on h and Rb + rc, not on the track radius.
%
%   This is the ONE place the mode -> surface/radius/sign mapping is
%   computed; every other file that needs it (rolling_matrices.m,
%   detect_events.m, handle_impact.m, lqr_design.m, linearize_system.m,
%   ball_hoop_dynamics.m, cascade_dynamics_reduced.m,
%   plan_trajectory_casadi.m) calls this function instead of
%   re-deriving it locally (docs/AUDIT.md item A3; formerly
%   mode_geometry.m, renamed and consolidated further per Phase 2).

    % Vee geometry, shared by all three tracks: tangency between the ball
    % and both O-ring tori places the ball centre a radial distance Delta
    % from the rail centreline, and puts the two contact points on a circle
    % of radius r_roll about the ball's spin axis.
    ball_radius = params.ball_radius;
    cord_radius = params.oring_cord_radius;
    half_gap    = params.oring_axial_spacing / 2;

    if half_gap >= ball_radius + cord_radius
        error('hoop_geometry:rails_too_far', ...
            ['O-ring separation %.2f mm exceeds 2*(ball_radius + cord_radius) = %.2f mm: ' ...
             'the ball would fall between the rails.'], ...
            2e3 * half_gap, 2e3 * (ball_radius + cord_radius));
    end

    Delta  = sqrt((ball_radius + cord_radius)^2 - half_gap^2);
    r_roll = ball_radius * Delta / (ball_radius + cord_radius);

    switch mode
        case 'rolling_out'
            surface = 'outer';
            R_track = params.outer_hoop_oring_radius;
            concave = true;
            psi_eq  = 0;

        case 'rolling_in_outside'
            % The inner hoop's rails sit on its rod surfaces, so its
            % outer/inner radii already ARE the rail centreline radii --
            % unlike the outer hoop, whose rails stand proud of its
            % rolling surface and need their own measured radius.
            surface = 'inner_outside';
            R_track = params.inner_hoop_outer_radius;
            concave = false;
            psi_eq  = pi;

        case 'rolling_in_inside'
            surface = 'inner_inside';
            R_track = params.inner_hoop_inner_radius;
            concave = true;
            psi_eq  = 0;

        otherwise
            error('hoop_geometry: unknown mode ''%s''.', mode);
    end

    if concave
        coupling_sign = -1;
        R_eff         = R_track - Delta;
    else
        coupling_sign = +1;
        R_eff         = R_track + Delta;
    end
end
