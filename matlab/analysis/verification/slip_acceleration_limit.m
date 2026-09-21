function [u_slip, info] = slip_acceleration_limit(params, mode)
%SLIP_ACCELERATION_LIMIT  Largest hoop acceleration the ball can follow
%                         without slipping, at the bottom of the track,
%                         at rest.
%
%   [u_slip, info] = slip_acceleration_limit(params, mode)
%
%   Outputs
%   -------
%   u_slip : bound on |theta_ddot|                              [rad/s^2]
%   info   : struct with .k (effective inertia factor), .R_track,
%            .mu, and .u_torque_equivalent [N.m] -- the torque that
%            acceleration corresponds to, for comparison against
%            cfg.tau_max
%
%   Derivation. A sphere carried by a surface whose contact point
%   accelerates at a_s needs a friction force F = m*a_s*(1 - 1/k) to roll
%   with it, where k = 1 + I_ball/(m*r_roll^2) is the effective inertia
%   factor of the two-rail contact (hoop_geometry.m; k = 1.552 here, not
%   the 1.4 of a single-point contact). With a_s = R_track*theta_ddot and
%   the available friction capped at mu*N,
%
%       |theta_ddot| <= mu*N / (m * (1 - 1/k) * R_track),
%
%   and at the bottom at rest N = m*g, giving the constant this function
%   returns.
%
%   THIS IS A SIZING NUMBER, NOT A TRAJECTORY CONSTRAINT. It holds only
%   where N = m*g. Along any real T3/T4 trajectory N varies by more than
%   an order of magnitude -- larger at the bottom of a fast swing
%   (centripetal term), and falling to zero at the top, where the
%   admissible acceleration falls to zero with it. Use
%   dynamics/reset_maps/contact_forces.m pointwise for anything that must actually
%   hold along a trajectory; use this one to size a drive amplitude or to
%   quote a single figure.
%
%   mu is params.ball_track_friction_coeff, which is ASSUMED and not
%   measured -- see the TODO(measure) note next to it. Anything reported
%   from this function should be reported with a sweep over mu, not at the
%   nominal value alone.

    [~, ~, R_track, ~, ~, r_roll] = hoop_geometry(mode, params);

    k  = 1 + params.ball_inertia / (params.ball_mass * r_roll^2);
    mu = params.ball_track_friction_coeff;

    u_slip = mu * params.gravity / ((1 - 1/k) * R_track);

    info.k                     = k;
    info.R_track               = R_track;
    info.mu                    = mu;
    info.u_torque_equivalent   = u_slip * params.total_inertia;
end
