function [A, B, C, D] = linearize_system(psi_eq, params, mode)
%LINEARIZE_SYSTEM  Linearize the ball-hoop dynamics around (psi_eq, 0)
%                  in CASCADE formulation (u = theta_ddot as input).
%
%   [A, B, C, D] = linearize_system(psi_eq, params, mode)
%
%   Produces the state-space matrices of the linearized system, valid
%   for small perturbations around the equilibrium (q_eq, q_dot_eq = 0).
%   The input is the angular acceleration u = theta_ddot, as in
%   Gurtner & Zemanek (2017).
%
%   State:  dx = [d_theta; d_theta_dot; d_psi; d_psi_dot]
%   Input:  du = d_theta_ddot                              [rad/s^2]
%   Output: dy = d_psi                                     [rad]
%
%   Inputs:
%     psi_eq : equilibrium ball position [rad]
%               - 'rolling_out'        : use 0 (stable)
%               - 'rolling_in_outside' : use pi (inverted, stabilizable)
%               - 'rolling_in_inside'  : use 0 (stable)
%     params : physical parameters from ball_hoop_params()
%     mode   : 'rolling_out' | 'rolling_in_outside' | 'rolling_in_inside'
%
%   Outputs:
%     A : 4x4 state matrix
%     B : 4x1 input matrix
%     C : 1x4 output matrix (selecting psi)
%     D : 1x1 feedthrough (zero)
%
%   Theory:
%     Starting from the full Lagrangian equation
%       M * [theta_ddot; psi_ddot] + C * [theta_dot; psi_dot] + G = [tau; 0]
%     in cascade mode we discard the theta-equation (handled by the inner
%     motor loop) and substitute theta_ddot = u in the psi-equation:
%       M(2,1)*u + M(2,2)*psi_ddot + C(2,1)*theta_dot + C(2,2)*psi_dot
%                                  + m*g*R_eff*sin(psi) = 0
%     Linearizing around (psi_eq, 0):
%       d_psi_ddot = (-M(2,1)*du - C(2,1)*d_theta_dot - C(2,2)*d_psi_dot
%                                - m*g*R_eff*cos(psi_eq)*d_psi) / M(2,2)
%
%   References:
%     Gurtner & Zemanek (2017), "Ball in double hoop", IFAC, eq. 7.

    %% --- Hoop geometry for this mode ---
    [~, R_eff] = hoop_geometry(mode, params);

    %% --- Evaluate dynamic matrices at the equilibrium ---
    % rolling_matrices already returns the 2x2 reduced (theta, psi) form
    % with the radial DOF eliminated by the rolling constraint.
    [M_eq, C_eq, ~] = rolling_matrices(params, psi_eq, mode);
    
    %% --- Gravity stiffness term ---
    % Linearization of m * g * R_eff * sin(psi) around psi_eq gives the
    % stiffness gravity_stiffness = m * g * R_eff * cos(psi_eq).
    %   - cos(0)  = +1  -> positive stiffness  (stable equilibrium)
    %   - cos(pi) = -1  -> negative stiffness  (unstable equilibrium)
    gravity_stiffness = params.ball_mass * params.gravity * R_eff * cos(psi_eq);

    %% --- Assemble state-space matrices ---
    % Extract relevant entries from the reduced 2x2 system.
    ball_effective_inertia      = M_eq(2,2);
    ball_hoop_inertial_coupling = M_eq(2,1);
    ball_hoop_damping_coupling  = C_eq(2,1);
    ball_damping                = C_eq(2,2);

    % Row 4 of the linearized ball equation, divided through by the ball's
    % own effective inertia (see the "Theory" section above).
    d_psi_ddot_d_theta_dot = -ball_hoop_damping_coupling / ball_effective_inertia;
    d_psi_ddot_d_psi       = -gravity_stiffness           / ball_effective_inertia;
    d_psi_ddot_d_psi_dot   = -ball_damping                / ball_effective_inertia;
    d_psi_ddot_d_u         = -ball_hoop_inertial_coupling  / ball_effective_inertia;

    % Row 1: d/dt(theta)     = theta_dot
    % Row 2: d/dt(theta_dot) = u                    (cascade: B(2,1) = 1)
    % Row 3: d/dt(psi)       = psi_dot
    % Row 4: d/dt(psi_dot)   = (linearized ball eq, with u substituted)

    A = [ 0, 1,                       0,                 0                     ;
          0, 0,                       0,                 0                     ;
          0, 0,                       0,                 1                     ;
          0, d_psi_ddot_d_theta_dot,  d_psi_ddot_d_psi,  d_psi_ddot_d_psi_dot ];

    B = [ 0             ;
          1             ;
          0             ;
          d_psi_ddot_d_u ];
    
    % Output: we measure psi (the 3rd state).
    C = [0, 0, 1, 0];
    D = 0;
end