function [new_mode, x0] = handle_impact(x_pre, params, cfg)
%HANDLE_IMPACT  Applies no-slip impact dynamics when the ball contacts a hoop
%               during free fall.
%
%   [new_mode, x0] = handle_impact(x_pre, params, cfg)
%
%   Inputs
%   ------
%   x_pre  : state vector just before impact [r; r_dot; theta; theta_dot; psi; psi_dot]
%   params : physical parameters struct (see ball_hoop_params)
%   cfg    : numerical settings struct (see sim_config)
%
%   Outputs
%   -------
%   new_mode : mode string after impact (or 'free_fall' if ball passes through gap)
%   x0       : initial state for the next integration segment

ball_r     = x_pre(1);
hoop_theta = x_pre(3);
ball_psi   = x_pre(5);

ball_psi_norm = mod(ball_psi, 2*pi);

% Hole boundaries in the inertial frame, at the current hoop rotation --
% must match detect_events.m's test exactly, since both decide the same
% physical question (is the ball over the gap) from the same event.
[hole_start, hole_end] = hole_bounds(params, hoop_theta);

% Cone of "free-fall". Wraparound-safe: the hole can straddle 0/2*pi once
% hoop_theta is nonzero, same as detect_events.m.
if hole_start < hole_end
    in_hole = (ball_psi_norm >= hole_start - cfg.TOL_ANGLE) && ...
              (ball_psi_norm <= hole_end   + cfg.TOL_ANGLE);
else
    in_hole = (ball_psi_norm >= hole_start - cfg.TOL_ANGLE) || ...
              (ball_psi_norm <= hole_end   + cfg.TOL_ANGLE);
end

% The gap belongs to the INNER hoop, so whether the ball is over it can
% only excuse an impact ON THE INNER HOOP. This test used to run here,
% before the radius was looked at, and returned free_fall for ANY contact
% whose azimuth happened to fall in the gap sector -- including a landing
% on the outer hoop, which is solid over all 360 deg. The ball then fell
% straight through it. Visible in the animation as the ball crossing the
% outer hoop instead of bouncing on it; it also let a T4 run bounce in
% and out of the inner hoop repeatedly, because each re-entry was
% adjudicated by the wrong test.
%
% The check is therefore applied inside the inner-hoop branch below, not
% before the branch. Nothing else about it changed.

% Determine contact surface by comparing radial position to the ball-centre
% orbit radius of each candidate mode. These must be hoop_geometry.m's
% values and not surface_radius -/+ ball_radius, both because the two-rail
% tracks put contact elsewhere and because detect_events.m fires on exactly
% these radii -- if the two disagreed by more than cfg.TOL_RADIUS the
% landing would be detected and then rejected, and the ball would fall
% straight through the hoop.
[~, R_eff_in_inside]  = hoop_geometry('rolling_in_inside',  params);
[~, R_eff_in_outside] = hoop_geometry('rolling_in_outside', params);
[~, R_eff_out]        = hoop_geometry('rolling_out',        params);

if abs(ball_r - R_eff_in_inside)  < cfg.TOL_RADIUS || ...
   abs(ball_r - R_eff_in_outside) < cfg.TOL_RADIUS
    % We're near the inner hoop, on one of its sides -- and this is the
    % one place where being over the gap means there is nothing to hit.
    if in_hole
        new_mode = 'free_fall';
        x0 = x_pre;
        return;
    end

    if ball_r < params.inner_hoop_inner_radius
        % We are inside the inner hoop
        new_mode  = 'rolling_in_inside';
    else
        new_mode  = 'rolling_in_outside';
    end
    [~, R_eff, R_track, coupling_sign, ~, r_roll] = hoop_geometry(new_mode, params);

    x0    = x_pre;
    x0(1) = R_eff;                                              % snap to surface [m]
    x0(2) = 0;                                                  % zero radial velocity [m/s]
    x0(5) = ball_psi_norm;                                      % normalise angle [rad]
    x0(6) = post_impact_psi_dot(x_pre(6), x0(4), R_eff, R_track, r_roll, coupling_sign, params);

elseif abs(ball_r - R_eff_out) < cfg.TOL_RADIUS
    % We're near the outer hoop, on its inner side only

    new_mode = 'rolling_out';
    [~, R_eff, R_track, coupling_sign, ~, r_roll] = hoop_geometry(new_mode, params);

    x0    = x_pre;
    x0(1) = R_eff;
    x0(2) = 0;
    x0(5) = ball_psi_norm;
    x0(6) = post_impact_psi_dot(x_pre(6), x0(4), R_eff, R_track, r_roll, coupling_sign, params);

else
    new_mode = 'free_fall';
    x0 = x_pre;
end
end


function psi_dot_post = post_impact_psi_dot(psi_dot_pre, theta_dot, R_eff, R_track, r_roll, coupling_sign, params)
%POST_IMPACT_PSI_DOT  Orbital rate of the ball just after it lands and
%                     starts rolling, from angular momentum about the
%                     contact line.
%
%   psi_dot_post = post_impact_psi_dot(psi_dot_pre, theta_dot, R_eff, R_track, r_roll, params)
%
%   THE LAW
%       psi_dot+ = (1/k)*psi_dot-
%                  + ((k-1)/k) * (r_roll*phi_dot- - cs*R_track*theta_dot) / R_eff
%       k  = 1 + I_ball/(m*r_roll^2)         (the two-rail factor, 1.5522)
%       cs = coupling_sign, -1 concave, +1 convex
%   with phi_dot- the ball's own spin, carried across the flight by
%   ball_spin_store.m.
%
%   DERIVATION. The impact is impulsive, so only the impulses matter.
%   Both act at the contact: the friction impulse along the contact line
%   and the normal impulse on a line that meets it. Neither therefore has
%   a moment about the CONTACT LINE, so the ball's angular momentum about
%   that line is conserved through the impact. Writing v = R_eff*psi_dot
%   for the ball-centre tangential speed and imposing rolling afterwards,
%
%       I_b*r_roll*phi_dot-  +  m*r_roll^2*v-  +  I_b*u_eff
%           = (I_b + m*r_roll^2) * v+ ,        u_eff = -cs*R_track*theta_dot
%
%   and dividing by (I_b + m*r_roll^2) = m*r_roll^2*k gives the law above.
%   It is a weighted average of what the ball brings and what the track
%   imposes, the weights being its translational against its rotational
%   inertia, so it is a convex combination and cannot run away.
%
%   WHY u_eff CARRIES cs, AND HOW THAT SIGN IS PINNED. The rolling
%   relation this codebase already uses (ball_spin_from_rolling.m,
%   contact_forces.m) is phi_dot = cs*(R_track/r_roll)*theta_dot
%   + (R_eff/r_roll)*psi_dot, which rearranges to
%   v - phi_dot*r_roll = -cs*R_track*theta_dot. So the track speed the
%   contact must match is -cs*R_track*theta_dot, not +R_track*theta_dot.
%   The sign is therefore not chosen here; it is inherited from the
%   rolling constraint. It only became visible once phi_dot- stopped
%   being zero -- with phi_dot- = 0 the two conventions differ by a term
%   that vanishes, which is how the earlier version of this function got
%   away with omitting cs entirely on the convex track.
%
%   THREE CHECKS IT PASSES. (i) A ball already rolling with the track is
%   left exactly unchanged -- substitute phi_dot- from the rolling
%   relation and the bracket collapses to R_eff*psi_dot-, giving
%   psi_dot+ = psi_dot-. (ii) On a STATIONARY track with no spin it keeps
%   1/k of its speed, E+/E- = 1/k = 0.644, the textbook sliding-to-rolling
%   loss, dissipative and independent of the incoming speed. (iii) Both
%   are checked numerically in analysis/studies/T4/t4c_impact_check.m.
%
%   WHAT IT REPLACES, AND WHY THAT WAS WRONG.
%       psi_dot+ = coupling_sign * (R_eff/r_roll) * theta_dot
%   took the post-impact ORBITAL rate from the hoop's speed alone and
%   discarded the ball's incoming tangential velocity entirely. With the
%   hoop at rest it set psi_dot+ = 0, so a ball landing tangentially on a
%   stationary track stopped dead whatever its speed -- E+/E- = 0, which
%   no impact law can produce. Its shape is the no-slip SPIN relation of
%   contact_forces.m written into the psi_dot slot: it is the ball's spin
%   that the track speed sets on contact, not its orbital rate. The gain
%   it applied to theta_dot, R_eff/r_roll = 3.06 on the inner-inside
%   track, is 6.3 times the correct ((k-1)/k)*(R_track/R_eff) = 0.486.
%   analysis/studies/T4/t4c_impact_check.m is the test.
%
%   THE SPIN IS NOT ASSUMED ZERO ANY MORE. It is locked to the other
%   states while the ball rolls, but frozen and independent once
%   airborne, so it cannot be recovered from the state AT IMPACT -- only
%   from the state at LIFTOFF, where the rolling relation still holds.
%   transition_state.m computes it there and ball_spin_store.m carries
%   it. Nothing integrates it: in flight no torque acts about the ball's
%   centre (gravity acts at the centre; drag is 1e-4 of the weight), so
%   it is a constant to remember, not a seventh state to propagate.
%   It is NOT a small term: on the T4c release at 125 deg the liftoff
%   spin is -192.6 rad/s and I_b*r_roll*phi_dot- is the same order as
%   m*r_roll^2*v-.
%
%   REMAINING ASSUMPTION -- the hoop is infinitely massive. The tangential
%   impulse acts on the hoop too, so theta_dot should drop by
%   J*R_track/J_hoop; this function does not change it. The hoop carries
%   J_hoop = 1.12e-3 kg.m^2 against the ball's m*R_eff^2 = 7.1e-5 on the
%   inner-inside track, a ratio of ~16, so the approximation is
%   reasonable there and worse on the outer track (ratio ~2). Beyond
%   ~10 % accuracy on the post-impact state, the coupled three-equation
%   impact (ball about the contact line, total angular momentum about the
%   hoop axis, rolling afterwards) is the version to write.

    k = 1 + params.ball_inertia / (params.ball_mass * r_roll^2);

    % The spin the ball carried into the flight, frozen since liftoff.
    phi_dot_pre = ball_spin_store('get');

    psi_dot_post = psi_dot_pre / k ...
                 + ((k-1)/k) * (r_roll * phi_dot_pre ...
                                - coupling_sign * R_track * theta_dot) / R_eff;
end
