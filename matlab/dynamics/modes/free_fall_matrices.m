function [M, C, G] = free_fall_matrices(r, psi, psi_dot, ...
        ball_mass, hoop_motor_inertia, motor_friction, gravity)
%FREE_FALL_MATRICES  Mass, Coriolis/damping, and gravity matrices for
%                    ball in free motion (no contact with either hoop).
%
%   [M, C, G] = free_fall_matrices(r, psi, psi_dot,
%       ball_mass, hoop_motor_inertia, motor_friction, gravity)
%
%   Inputs
%   ------
%   r                  : ball radial position [m]
%   psi                : ball angular position, inertial frame, from
%                        the downward vertical [rad]
%   psi_dot            : ball angular velocity [rad/s]
%   ball_mass          : ball mass [kg]
%   hoop_motor_inertia : combined hoop+motor rotational inertia [kg.m^2]
%   motor_friction     : motor viscous friction coefficient [N.m.s/rad]
%   gravity            : gravitational acceleration [m/s^2]
%
%   Outputs
%   -------
%   M : 3x3 mass/inertia matrix for the (r, theta, psi) system [kg, kg.m^2]
%   C : 3x3 Coriolis/damping matrix  (state-dependent)
%   G : 3x1 gravity vector                                     [N, N.m]
%
%   In free fall the radial DOF is unconstrained, giving a 3-DOF system.
%   The C matrix is velocity-dependent because of the centrifugal and
%   Coriolis terms from the polar-coordinate Lagrangian.

M = diag([ball_mass, hoop_motor_inertia, ball_mass * r^2]);

C = zeros(3, 3);
C(1, 3) = -ball_mass * r * psi_dot;
C(2, 2) =  motor_friction;
C(3, 1) =  2 * ball_mass * r * psi_dot;

G = [-ball_mass * gravity * cos(psi);    % radial gravity component
      0;                                 % hoop is torque-driven, no gravity term
      ball_mass * gravity * r * sin(psi)];
end
