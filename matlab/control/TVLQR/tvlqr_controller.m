function [u, dz_dt, info] = tvlqr_controller(t, x_phys, ref_traj, ctrl_params)
%TVLQR_CONTROLLER  Time-varying LQR with feedforward, applied along a
%                  reference trajectory.
%
%   [u, dz_dt, info] = tvlqr_controller(t, x_phys, ref_traj, ctrl_params)
%
%   Control law: u(t) = u*(t) - K(t) * (x(t) - x*(t))
%
%   Inputs
%   ------
%   t           : current time [s]
%   x_phys      : physical state [r; r_dot; theta; theta_dot; psi; psi_dot]
%                 [m; m/s; rad; rad/s; rad; rad/s]
%   ref_traj    : struct with fields:
%                   .t_traj  : time grid (N x 1)                     [s]
%                   .x_traj  : reference reduced state (N x 4)
%                              [theta, theta_dot, psi, psi_dot], [rad, rad/s, rad, rad/s]
%                   .u_traj  : reference control (N x 1)             [rad/s^2]
%                   .K_traj  : time-varying gain (N x 4)             [(rad/s^2)/state-unit, per column]
%                   .Tf      : final time of the trajectory          [s]
%   ctrl_params : struct with fields
%                   .tau_max        : torque saturation limit    [N.m]
%                   .total_inertia  : motor+hoop inertia         [kg.m^2]
%                   .motor_friction : viscous friction coeff.    [N.m.s/rad]
%
%   Outputs
%   -------
%   u      : commanded acceleration [rad/s^2]
%   dz_dt  : zeros(0,1) (no internal state)
%   info   : diagnostics

    %% --- After trajectory ends, hand off to a stationary regulator ---
    % Once t > Tf, hold the final state target using K_traj(end,:).
    if t <= ref_traj.Tf
        x_ref = interp1(ref_traj.t_traj, ref_traj.x_traj, t, 'linear', 'extrap').';
        u_ref = interp1(ref_traj.t_traj, ref_traj.u_traj, t, 'linear', 'extrap');
        K     = interp1(ref_traj.t_traj, ref_traj.K_traj, t, 'linear', 'extrap');
    else
        x_ref = ref_traj.x_traj(end, :).';
        u_ref = ref_traj.u_traj(end);
        K     = ref_traj.K_traj(end, :);
    end
    
    %% --- Extract reduced state from physical state ---
    % Physical state: [r; r_dot; theta; theta_dot; psi; psi_dot]
    % Reduced state: [theta; theta_dot; psi; psi_dot]
    x_red = [x_phys(3); x_phys(4); x_phys(5); x_phys(6)];
    
    %% --- TVLQR control law with feedforward ---
    dx = x_red - x_ref;
    u_unsat = u_ref - K * dx;
       
    %% --- Saturations physiques ---------------------------------------------

    theta_dot = x_phys(4);
    
    % Limite d'acceleration due au couple moteur disponible.
    u_motor_max = ...
        ( ctrl_params.tau_max ...
        - ctrl_params.motor_friction * theta_dot) ...
        / ctrl_params.total_inertia;
    
    u_motor_min = ...
        (-ctrl_params.tau_max ...
        - ctrl_params.motor_friction * theta_dot) ...
        / ctrl_params.total_inertia;
    
    % Limite supplementaire de non-glissement.
    u_hi = min(u_motor_max,  ctrl_params.u_slip_max);
    u_lo = max(u_motor_min, -ctrl_params.u_slip_max);
    
    u = max(u_lo, min(u_hi, u_unsat));

    %% --- Protection de la limite de vitesse -------------------------------

    theta_dot_max = ctrl_params.theta_dot_max;
    
    if theta_dot >= theta_dot_max && u > 0
        u = 0;
    elseif theta_dot <= -theta_dot_max && u < 0
        u = 0;
    end

    %% --- No internal state ---
    dz_dt = zeros(0, 1);
    
    %% --- Diagnostics ---
    info.x_ref     = x_ref;
    info.u_ref     = u_ref;
    info.dx        = dx;
    info.u_unsat = u_unsat;
    info.u       = u;
    info.theta_dot = theta_dot;
    info.u_motor_max = u_motor_max;
    info.u_motor_min = u_motor_min;
    info.saturated = abs(u - u_unsat) > 1e-12;
    info.speed_limited = ...
        (theta_dot >= theta_dot_max && u_unsat > 0) || ...
        (theta_dot <= -theta_dot_max && u_unsat < 0);
end