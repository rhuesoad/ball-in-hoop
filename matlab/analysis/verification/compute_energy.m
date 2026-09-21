function [E_total, E_kin, E_pot, E_error_percent] = compute_energy(t, X, params, mode_history)
%COMPUTE_ENERGY  Mechanical energy of the ball-hoop system along a trajectory.
%
%   [E_total, E_kin, E_pot, E_error_percent] = compute_energy(t, X, params, mode_history)
%
%   Inputs
%   ------
%   t            : time vector [s]  (N x 1)
%   X            : state matrix (N x 6)
%                  columns: [r, r_dot, theta, theta_dot, psi, psi_dot]
%   params       : physical parameters struct (see ball_hoop_params)
%   mode_history : (optional) mode string at each time step — reserved for
%                  future use; currently unused
%
%   Outputs
%   -------
%   E_total         : total mechanical energy at each step [J]  (N x 1)
%   E_kin           : kinetic energy [J]  (N x 1)
%   E_pot           : gravitational potential energy [J]  (N x 1)
%   E_error_percent : drift relative to E_total(1) [%]  (N x 1)
%                     Non-zero because of motor energy input and numerical
%                     dissipation; use as a solver quality indicator.

%#ok<INUSD> mode_history reserved for per-mode spin velocity formulas

ball_r         = X(:, 1);
ball_r_dot     = X(:, 2);
ball_psi       = X(:, 5);
ball_psi_dot   = X(:, 6);
hoop_theta_dot = X(:, 4);

% Radius of the ball centre when in rolling_out contact, and the radius it
% actually rolls on (two-rail O-ring track -- see hoop_geometry.m). r_roll
% is smaller than the ball radius, so the ball spins faster than the naive
% no-slip formula predicts and stores more energy in spin; using
% ball_radius here understates T_ball_rot by about 38 %.
[~, R_eff_out, R_track_out, ~, ~, r_roll] = hoop_geometry('rolling_out', params);

% Tolerance for detecting rolling contact: ball centre within 1 mm of R_eff.
% Used only in the spin-velocity formula; does not affect E_kin or E_pot.
ROLLING_CONTACT_TOL = 1e-3;                         % [m]

n     = length(t);
E_kin = zeros(n, 1);
E_pot = zeros(n, 1);

for i = 1:n
    %% Kinetic energy — hoop
    T_hoop = 0.5 * params.hoop_motor_inertia * hoop_theta_dot(i)^2;

    %% Kinetic energy — ball (translational)
    v_ball_sq    = ball_r_dot(i)^2 + (ball_r(i) * ball_psi_dot(i))^2;
    T_ball_trans = 0.5 * params.ball_mass * v_ball_sq;

    %% Kinetic energy — ball (spin)
    if abs(ball_r(i) - R_eff_out) < ROLLING_CONTACT_TOL
        % Rolling contact: no-slip constraint gives the ball spin rate.
        % phi_dot = (R_track / r_roll) * theta_dot - (R_eff / r_roll) * psi_dot
        radius_ratio = R_track_out / r_roll;
        phi_dot      = radius_ratio * hoop_theta_dot(i) ...
                     - (R_eff_out / r_roll) * ball_psi_dot(i);
        T_ball_rot = 0.5 * params.ball_inertia * phi_dot^2;
    else
        % Free-fall or inner-hoop contact: spin kinetic energy not computed
        % (no-slip formula not valid; spin is not tracked in the state vector).
        T_ball_rot = 0;
    end

    E_kin(i) = T_hoop + T_ball_trans + T_ball_rot;

    %% Potential energy
    % Reference: hoop centre. Negative at the bottom (psi=0) where r*cos(psi)>0.
    E_pot(i) = -params.ball_mass * params.gravity * ball_r(i) * cos(ball_psi(i));
end

E_total = E_kin + E_pot;

E_0 = E_total(1);
if abs(E_0) > 1e-10
    E_error_percent = (E_total - E_0) / abs(E_0) * 100;
else
    E_error_percent = zeros(n, 1);
end
end
