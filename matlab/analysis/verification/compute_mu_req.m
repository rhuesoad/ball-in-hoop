function [mu_req, N, F_t] = compute_mu_req(x_traj, u_traj, params, mode)
%COMPUTE_MU_REQ  Friction coefficient the contact must supply, node by node.
%
%   mu_req = compute_mu_req(x_traj, u_traj, params, mode)
%   [mu_req, N, F_t] = compute_mu_req(x_traj, u_traj, params, mode)
%
%   Inputs
%   ------
%   x_traj : reduced state trajectory (K x 4)
%            [theta, theta_dot, psi, psi_dot], [rad, rad/s, rad, rad/s]
%   u_traj : commanded hoop acceleration at the same nodes (K x 1) [rad/s^2]
%   params : physical parameters (see ball_hoop_params)
%   mode   : 'rolling_out' | 'rolling_in_outside' | 'rolling_in_inside'
%
%   Outputs
%   -------
%   mu_req : |F_t| / N at every node (K x 1)                           [-]
%   N      : normal contact force at every node (K x 1)                [N]
%   F_t    : tangential contact force at every node (K x 1)            [N]
%
%   WHAT THIS IS. mu_req is the smallest static friction coefficient for
%   which the plan is a ROLLING plan. Where mu_req exceeds the coefficient
%   the contact can actually supply, the trajectory is one the ball cannot
%   follow without slipping, and the whole reduced model -- which assumes
%   rolling -- stops describing it. It is the diagnostic form of the
%   constraint plan_trajectory_casadi.m can impose (opts.enforce_no_slip):
%   the constraint asks "is mu*N - |F_t| >= 0 for THIS mu", this function
%   asks "which mu would be needed", which is the question to ask when mu
%   itself has not been measured (ball_hoop_params.m marks
%   ball_track_friction_coeff CONSERVATIVE, not measured).
%
%   WHERE EACH TERM COMES FROM. Both forces are taken from
%   dynamics/reset_maps/contact_forces.m, which is the single derivation of them in
%   this codebase and is checked numerically against row 2 of
%   rolling_matrices.m in analysis/verification/run_all_checks.m (check9). Nothing is
%   re-derived here; this function only forms the ratio and handles the
%   nodes where the ratio is not defined. For the record, the two
%   equations contact_forces.m evaluates are:
%
%     N     radial equilibrium of the ball on its circular orbit,
%           N = -coupling_sign * m * (g*cos(psi) + R_eff*psi_dot^2),
%           the same expression detect_events.m uses for liftoff, with
%           coupling_sign and R_eff from hoop_geometry.m.
%
%     F_t   the ball's angular equation about its own centre,
%           I_ball*phi_ddot = -F_t*r_roll - b_ball*phi_dot,
%           with phi_ddot supplied by the derivative of the no-slip
%           constraint, phi_ddot = coupling_sign*rho*u
%           + (R_eff/r_roll)*psi_ddot, rho = R_track/r_roll, and psi_ddot
%           from cascade_dynamics_reduced.m.
%
%           NOTE THE LEVER ARM. It is r_roll, NOT the ball radius. The
%           ball rides in a vee formed by two O-rings and touches each
%           track at two points, so it spins about a chord of itself and
%           the contact points sit at r_roll = Rb*Delta/(Rb + rc) < Rb
%           from that axis (hoop_geometry.m). Writing I_ball*phi_ddot =
%           -F_t*Rb would be the single-point-contact form and is wrong
%           for this bench. The rolling-contact dissipation b_ball*phi_dot
%           acts on the ball's spin, so it belongs to this equation too.
%
%   WHERE THE RATIO IS NOT DEFINED. mu_req = |F_t|/N presumes a contact.
%   Near the top of a loop N collapses toward zero, and no finite friction
%   coefficient delivers a finite F_t there -- so the ratio is returned as
%   Inf for N <= 0, not clipped and not skipped. A margin sweep is exactly
%   the study in which those nodes appear, and reporting them as a large
%   finite number would hide the fact that the binding failure at zero
%   margin is loss of contact, not loss of grip.

    if size(x_traj, 2) ~= 4
        error('compute_mu_req:x_traj_shape', ...
            'x_traj must be K x 4 ([theta, theta_dot, psi, psi_dot]), got %d x %d.', ...
            size(x_traj, 1), size(x_traj, 2));
    end
    u_traj = u_traj(:);
    K = size(x_traj, 1);
    if numel(u_traj) ~= K
        error('compute_mu_req:length_mismatch', ...
            'u_traj has %d entries but x_traj has %d nodes.', numel(u_traj), K);
    end

    mu_req = zeros(K, 1);
    N      = zeros(K, 1);
    F_t    = zeros(K, 1);

    for i = 1:K
        [N(i), F_t(i)] = contact_forces(x_traj(i,:).', u_traj(i), params, mode);
        if N(i) > 0
            mu_req(i) = abs(F_t(i)) / N(i);
        else
            % No contact (or exactly on the verge): the ball is not being
            % held by friction at all, so no coefficient makes this node a
            % rolling node. See the docstring.
            mu_req(i) = Inf;
        end
    end
end
