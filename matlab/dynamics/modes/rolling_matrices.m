function [M, C, G] = rolling_matrices(params, psi, mode)
%ROLLING_MATRICES  Mass, Coriolis/damping, and gravity matrices for
%                  ball rolling on a hoop surface.
%
%   [M, C, G] = rolling_matrices(params, psi, mode)
%
%   Inputs
%   ------
%   params : physical parameters struct (see ball_hoop_params)
%   psi    : ball angular position, inertial frame                  [rad]
%   mode   : 'rolling_out' | 'rolling_in_outside' | 'rolling_in_inside'
%
%   Outputs
%   -------
%   M : 2x2 mass/inertia matrix for the (theta, psi) subsystem     [kg.m^2]
%   C : 2x2 damping matrix                                      [N.m.s/rad]
%   G : 2x1 gravity vector                                            [N.m]
%
%   The radial degree of freedom is eliminated by the rolling constraint,
%   leaving a 2-DOF system in (theta, psi). Equation of motion:
%       M * [theta_ddot; psi_ddot] + C * [theta_dot; psi_dot] + G = [tau; 0]
%
%   Surface geometry (R_eff, contact radius, coupling sign) comes from
%   hoop_geometry.m, the single source of truth for the mode -> surface
%   mapping (docs/AUDIT.md item A3; this function used to keep its own
%   independent switch on a lower-level surface_type argument -- Phase 2
%   consolidated it into hoop_geometry.m like every other caller).

%% --- Unpack physical parameters ---
    ball_mass          = params.ball_mass;
    ball_inertia       = params.ball_inertia;
    hoop_motor_inertia = params.hoop_motor_inertia;
    ball_friction      = params.ball_friction;
    motor_friction     = params.motor_friction;
    gravity            = params.gravity;

%% --- Surface-specific geometry ---
    [~, R_eff, R_track, coupling_sign, ~, r_roll] = hoop_geometry(mode, params);

% Every no-slip ratio below divides by r_roll, the radius the ball actually
% rolls on, NOT params.ball_radius. The two differ because the ball rests in
% the vee of two O-rings (hoop_geometry.m): r_roll = 10.64 mm against a
% 12.5 mm ball. ball_inertia keeps ball_radius -- it is the ball's own
% inertia (2/5)m*Rb^2, a property of the sphere, not of how it is carried.
% Conflating the two is the easy mistake here: the result stays plausible
% and is wrong by ~11 % on every coupling term.
radius_ratio = R_track / r_roll;   % no-slip angular velocity ratio

% Transferred to the hoops center
ball_effective_inertia = ball_mass * R_eff^2 + ball_inertia * (R_eff / r_roll)^2;
inertia_coupling = ball_inertia * radius_ratio * (R_eff / r_roll);     % From no-slip
M = [hoop_motor_inertia + ball_inertia * radius_ratio^2,  coupling_sign * inertia_coupling;
     coupling_sign * inertia_coupling,                    ball_effective_inertia];

% Dissipation
C = [motor_friction + ball_friction * radius_ratio^2, ...
     coupling_sign * ball_friction * radius_ratio * (R_eff / r_roll); ...
     coupling_sign * ball_friction * radius_ratio * (R_eff / r_roll), ...
     ball_friction * (R_eff / r_roll)^2];

% Gravity only on psi
G = [0; ball_mass * gravity * R_eff * sin(psi)];
end
