function [value, isterminal, direction] = detect_events(~, x, params, cfg, current_mode)
%DETECT_EVENTS  Zero-crossing event functions for ball-hoop mode transitions.
%
%   [value, isterminal, direction] = detect_events(t, x, params, cfg, current_mode)
%
%   Inputs
%   ------
%   t            : current time [s]  (unused, required by odeset Events API)
%   x            : state vector [r; r_dot; theta; theta_dot; psi; psi_dot]
%   params       : physical parameters struct (see ball_hoop_params)
%   cfg          : numerical settings struct (see sim_config)
%   current_mode : active dynamic mode string
%
%   Outputs (odeset Events convention)
%   -----------------------------------
%   value      : event function values — integration stops when any crosses zero
%   isterminal : 1 = stop integration at this event
%   direction  : 0 = detect any crossing direction

ball_r       = x(1);
hoop_theta   = x(3);
ball_psi     = x(5);
ball_psi_dot = x(6);

% Disable events if we are near the origin, because degenerate value here.
% Event-function output size must match the count each mode's branch
% below declares (odeset requires a fixed count within one integration
% segment) -- docs/AUDIT.md item A5.
switch current_mode
    case {'rolling_out', 'rolling_in_outside'}
        n_events = 1;
    case {'free_fall', 'rolling_in_inside'}
        n_events = 2;
    otherwise
        error('detect_events: unknown mode ''%s''.', current_mode);
end

if ball_r < params.ball_radius
    value      = ones(n_events, 1);
    isterminal = ones(n_events, 1);
    direction  = zeros(n_events, 1);
    return;
end

% Hole boundaries updated with movement of hoops
[hole_start, hole_end] = hole_bounds(params, hoop_theta);
hole_start = mod(hole_start - cfg.ANGLE_PAD, 2*pi);
hole_end   = mod(hole_end   + cfg.ANGLE_PAD, 2*pi);

ball_psi_norm = mod(ball_psi, 2*pi);

% Detection of the hole region
if hole_start < hole_end
    in_hole = (ball_psi_norm >= hole_start) && (ball_psi_norm <= hole_end);
else
    in_hole = (ball_psi_norm >= hole_start) || (ball_psi_norm <= hole_end);
end


switch current_mode
    case 'rolling_out'
        % Normal force at the outer hoop. Liftoff occurs when N drops to zero.
        [~, R_eff] = hoop_geometry(current_mode, params);
        N     = params.ball_mass * (params.gravity * cos(ball_psi) + R_eff * ball_psi_dot^2);
        value      = N;
        isterminal = 1;
        direction  = 0;

    case 'free_fall'
        % Contact radii come from hoop_geometry.m, never from
        % surface_radius -/+ ball_radius: with the two-rail O-ring tracks
        % the ball centre settles a distance Delta from the rail
        % centreline, not one ball radius from a surface. Re-deriving them
        % here is what previously let this function and handle_impact.m
        % disagree with rolling_matrices.m about where contact happens.
        if ball_r < params.inner_hoop_radius
            [~, R_eff_inner] = hoop_geometry('rolling_in_inside', params);
        else
            [~, R_eff_inner] = hoop_geometry('rolling_in_outside', params);
        end
        val_inner = ball_r - R_eff_inner;
        if in_hole
            val_inner = 1;  % Cancel inner-hoop detection when ball is over the gap.
        end

        % Contact with outer hoop.
        [~, R_eff_outer] = hoop_geometry('rolling_out', params);
        val_outer = ball_r - R_eff_outer;

        value      = [val_inner; val_outer];
        isterminal = [1; 1];
        direction  = [0; 0];

    case 'rolling_in_outside'
        % Normal force at the outer (convex) surface of the inner hoop.
        % Convex contact flips the sign relative to the concave case
        % above: the surface pushes the ball outward, not inward, so the
        % centripetal and gravity-outward terms both flip (see
        % docs/AUDIT.md, item 1.2, for the radial free-body derivation).
        [~, R_eff] = hoop_geometry(current_mode, params);
        N     = - params.ball_mass * (params.gravity * cos(ball_psi) + R_eff * ball_psi_dot^2);
        value      = N;
        isterminal = 1;
        direction  = 0;

    case 'rolling_in_inside'
        %   1. Normal force N = 0  -> liftoff from inner surface
        %   2. Signed angular distance to hole reaches zero from above ->
        %      ball arrives at hole opening; direction = -1 ensures the
        %      event only fires as ball enters (not exits) the gap window.
        %      (replaces the value=0 shortcut which never created a zero-crossing)
        [~, R_eff] = hoop_geometry(current_mode, params);
        N        = params.ball_mass * (params.gravity * cos(ball_psi) + R_eff * ball_psi_dot^2);
        ang_dist = signed_angular_distance_to_hole(ball_psi, hoop_theta, params);
        value      = [N; ang_dist];
        isterminal = [1; 1];
        direction  = [0; -1];

    otherwise
        error('detect_events: unknown mode ''%s''.', current_mode);
end

end

function d = signed_angular_distance_to_hole(ball_psi, hoop_theta, params)
% Positive when ball is outside the hole angular window, zero at boundary,
% negative inside. direction = -1 triggers only on entry (positive -> negative).
    hole_center = mod(params.hole_center_angle + hoop_theta, 2*pi);
    ball_norm   = mod(ball_psi, 2*pi);
    delta       = mod(ball_norm - hole_center + pi, 2*pi) - pi;   % wrapped to (-pi, pi]
    d = abs(delta) - params.hole_angular_width / 2;
end
