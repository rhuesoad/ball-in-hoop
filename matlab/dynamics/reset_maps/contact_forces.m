function [N, F_t, slip_margin] = contact_forces(x_red, u, params, mode)
%CONTACT_FORCES  Normal and tangential contact force on the ball while it
%                rolls, and the margin left before it slips.
%
%   [N, F_t, slip_margin] = contact_forces(x_red, u, params, mode)
%
%   Inputs
%   ------
%   x_red  : reduced state [theta; theta_dot; psi; psi_dot]
%            [rad; rad/s; rad; rad/s]  (theta itself is unused -- it is a
%            cyclic coordinate, cascade_dynamics_reduced.m:31)
%   u      : commanded hoop acceleration theta_ddot          [rad/s^2]
%   params : physical parameters (see ball_hoop_params)
%   mode   : 'rolling_out' | 'rolling_in_outside' | 'rolling_in_inside'
%
%   Outputs
%   -------
%   N           : normal force pressing ball and track together      [N]
%                 (>= 0 while in contact; N = 0 is liftoff)
%   F_t         : tangential (friction) force at the contact         [N]
%   slip_margin : mu*N - |F_t|                                       [N]
%                 (>= 0 while the no-slip assumption holds; < 0 means the
%                 model's own rolling constraint is being violated)
%
%   WHY THIS FUNCTION EXISTS. The whole model rests on rolling without
%   slipping (docs/MODEL.md sec. 3), and nothing anywhere checked that the
%   friction needed to enforce it is actually available. For T1/T2 that
%   was harmless -- the torques involved are small. For T3/T4 it is not:
%   the manoeuvres ask for hoop accelerations of the same order as the
%   bound below, so a trajectory can be perfectly feasible for the motor
%   and still be one the ball cannot follow.
%
%   NORMAL FORCE. Radial equilibrium of the ball on its circular orbit:
%       N = -coupling_sign * m * (g*cos(psi) + R_eff*psi_dot^2)
%   the same expression detect_events.m:64 uses for liftoff detection, with
%   the concave/convex sign taken from hoop_geometry.m rather than a local
%   switch (concave: coupling_sign = -1, the wall pushes inward).
%
%   TANGENTIAL FORCE. Friction is the only tangential force the track can
%   apply, and it is what spins the ball up. The ball's angular equation
%   about its own centre carries both it and the rolling-contact
%   dissipation of docs/MODEL.md sec. 4.2 (which acts on phi_dot, i.e. on
%   the ball's spin, so it belongs to this equation and not to the
%   tangential one):
%       I_ball * phi_ddot = -F_t * r_roll - b_ball * phi_dot,
%   with the no-slip constraint of docs/MODEL.md eq. (1) and its
%   derivative supplying
%       phi_dot  = coupling_sign * rho * theta_dot + (R_eff/r_roll) * psi_dot,
%       phi_ddot = coupling_sign * rho * u         + (R_eff/r_roll) * psi_ddot,
%   and rho = R_track/r_roll. Substituting into the ball's tangential
%   equation m*R_eff*psi_ddot = F_t - m*g*sin(psi) reproduces row 2 of
%   rolling_matrices.m exactly, term for term, damping included -- which
%   is what fixes the sign conventions above, and is checked numerically
%   in analysis/verification/run_all_checks.m rather than asserted here.
%
%   THE SPECIAL CASE EVERYONE QUOTES. At the bottom, at rest
%   (psi = psi_dot = 0), the two expressions collapse to
%       |u| <= mu*g / ((1 - 1/k) * R_track),   k = 1 + I_ball/(m*r_roll^2),
%   which is the bound the bench code applies as a constant
%   (t1_config_v3.py, U_SLIP_MAX) and which
%   analysis/verification/slip_acceleration_limit.m returns. That constant is only
%   valid there. Near the top of the outer hoop N falls toward zero, so
%   the available friction does too and the admissible u collapses with
%   it -- which is precisely the part of a looping trajectory where the
%   planner is otherwise most tempted to spend torque. Use this function,
%   not the constant, anywhere the ball is not near the bottom.
%
%   CasADi-safe: only +, -, *, /, sin, cos and abs on the inputs, so the
%   trajectory planner can impose slip_margin >= 0 symbolically at every
%   collocation point.

    [~, R_eff, R_track, coupling_sign, ~, r_roll] = hoop_geometry(mode, params);

    psi     = x_red(3);
    psi_dot = x_red(4);

    m       = params.ball_mass;
    g       = params.gravity;
    I_ball  = params.ball_inertia;
    b_ball  = params.ball_friction;

    %% --- Normal force ---
    N = -coupling_sign * m * (g * cos(psi) + R_eff * psi_dot^2);

    %% --- Tangential (friction) force ---
    dx       = cascade_dynamics_reduced(x_red, u, params, mode);
    psi_ddot = dx(4);

    rho      = R_track / r_roll;
    phi_dot  = coupling_sign * rho * x_red(2) + (R_eff / r_roll) * psi_dot;
    phi_ddot = coupling_sign * rho * u        + (R_eff / r_roll) * psi_ddot;
    F_t      = -(I_ball * phi_ddot + b_ball * phi_dot) / r_roll;

    %% --- How much friction is left ---
    slip_margin = params.ball_track_friction_coeff * N - abs(F_t);
end
