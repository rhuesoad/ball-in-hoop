function [u, dz_dt, info] = lqr_controller(t, x_phys, ref, ctrl_params)
%LQR_CONTROLLER  Cascade LQR with feedforward and torque-saturation
%                translated to acceleration saturation.
%
%   [u, dz_dt, info] = lqr_controller(t, x_phys, ref, ctrl_params)
%
%   Inputs
%   ------
%   t           : current time [s]
%   x_phys      : physical state [r; r_dot; theta; theta_dot; psi; psi_dot]
%   ref         : reference struct with handles .psi(t) [rad], .psi_dot(t) [rad/s]
%   ctrl_params : struct from lqr_design with fields:
%                   .K              (1x4) LQR gain matrix [(rad/s^2)/state-unit, per column]
%                   .u_eq           feedforward acceleration [rad/s^2]
%                   .tau_max        : torque saturation limit    [N.m]
%                   .total_inertia  : motor+hoop inertia         [kg.m^2]
%                   .motor_friction : viscous friction coeff.    [N.m.s/rad]
%
%   Outputs
%   -------
%   u      : commanded angular acceleration u = theta_ddot_ref [rad/s^2]
%   dz_dt  : 0 (no internal state in basic LQR)
%   info   : diagnostic struct
%
%   Control law
%   -----------
%   The state deviation is:
%     dx = [theta - 0;
%           theta_dot - 0;
%           psi - psi_ref(t);
%           psi_dot - psi_dot_ref(t)]
%   The control is:
%     u_unsat = u_eq - K * dx
%   then saturated to respect the physical torque limit tau_max.

    %% --- Extract relevant states ---
    theta     = x_phys(3);
    theta_dot = x_phys(4);
    psi       = x_phys(5);
    psi_dot   = x_phys(6);
    
    %% --- Evaluate reference at current time ---
    psi_ref     = ref.psi(t);
    psi_dot_ref = ref.psi_dot(t);
    
    %% --- State deviation vector ---
    % theta and theta_dot reference is 0 (the hoop has no absolute reference
    % position; we only want to penalize unnecessary rotation).
    dx = [ theta;
           theta_dot;
           psi     - psi_ref;
           psi_dot - psi_dot_ref ];
    
    %% --- Feedforward acceleration ---
    % Static case (stabilization): u_eq compensates the gravity bias at the
    % operating point (constant, computed offline by lqr_design).
    % Dynamic case (tracking): if the reference provides psi_ddot, compute the
    % exact open-loop acceleration u_ff(t) that realizes the trajectory, by
    % inverting the ball equation INCLUDING the live Coriolis coupling
    % C(2,1)*theta_dot. The LQR then only corrects the deviation.
    if isfield(ref, 'psi_ddot') && ~isempty(ref.psi_ddot)
        psi_ddot_ref = ref.psi_ddot(t);
        u_ff = feedforward(psi_ref, psi_dot_ref, psi_ddot_ref, ...
                                 theta_dot, ctrl_params);
    else
        u_ff = ctrl_params.u_eq;
    end
    
    %% --- LQR law with feedforward ---
    u_unsat = u_ff - ctrl_params.K * dx;
    
    %% --- Translate torque saturation into acceleration saturation ---
    % tau = I_total * u + b_motor * theta_dot
    % => |tau| <= tau_max constrains u within bounds depending on theta_dot
    u_max = ( ctrl_params.tau_max - ctrl_params.motor_friction * theta_dot) ...
             / ctrl_params.total_inertia;
    u_min = (-ctrl_params.tau_max - ctrl_params.motor_friction * theta_dot) ...
             / ctrl_params.total_inertia;
    
    u = max(u_min, min(u_max, u_unsat));
    
    %% --- No internal state to integrate ---
    dz_dt = [];
    
    %% --- Diagnostics ---
    info.dx        = dx;
    info.u_unsat   = u_unsat;
    info.u         = u;
    info.psi_ref   = psi_ref;
    info.e_psi     = psi_ref - psi;
    info.saturated = (u ~= u_unsat);
end