function dx = ball_hoop_dynamics(t, x, tau_fun, params, mode)
%BALL_HOOP_DYNAMICS  State derivative for the ball-in-double-hoop system.
%
%   dx = ball_hoop_dynamics(t, x, tau_fun, params, mode)
%
%   Inputs
%   ------
%   t       : current time [s]
%   x       : state vector [r; r_dot; theta; theta_dot; psi; psi_dot]
%               r         = ball radial position [m]      (free-fall only)
%               r_dot     = ball radial velocity [m/s]    (free-fall only)
%               theta     = hoop angle [rad]
%               theta_dot = hoop angular velocity [rad/s]
%               psi       = ball angular position, inertial frame, from
%                           the downward vertical [rad]
%               psi_dot   = ball angular velocity [rad/s]
%   tau_fun : function handle  tau_fun(t, x) -> torque [N.m]
%   params  : physical parameters struct (see ball_hoop_params)
%   mode    : current dynamic mode — one of:
%               'free_fall'           ball in free motion
%               'rolling_out'         ball on inner surface of outer hoop
%               'rolling_in_outside'  ball on outer surface of inner hoop
%               'rolling_in_inside'   ball on inner surface of inner hoop
%
%   Outputs
%   -------
%   dx : state derivative [r_dot; r_ddot; theta_dot; theta_ddot; psi_dot; psi_ddot]

r       = x(1);   r_dot     = x(2);
theta   = x(3);   theta_dot = x(4);  
psi     = x(5);   psi_dot   = x(6);

tau = tau_fun(t, x);

switch mode
    case 'free_fall'
        [M, C, G] = free_fall_matrices(r, psi, psi_dot, ...
            params.ball_mass, params.hoop_motor_inertia, ...
            params.motor_friction, params.gravity);
        % 3-DOF system: r, theta, psi all free.
        generalized_velocities = [r_dot; theta_dot; psi_dot];
        B_input = [0; 1; 0];    % torque enters theta equation (index 2 of 3)

    case {'rolling_out', 'rolling_in_outside', 'rolling_in_inside'}
        [M, C, G] = rolling_matrices(params, psi, mode);
        % 2-DOF system: rolling constraint fixes r, leaving (theta, psi).
        generalized_velocities = [theta_dot; psi_dot];
        B_input = [1; 0];       % torque enters theta equation (index 1 of 2)

    otherwise
        error('ball_hoop_dynamics: unknown mode ''%s''.', mode);
end

% M * qddot = B * tau - C * qdot - G
generalized_accelerations = M \ (B_input * tau - C * generalized_velocities - G);

if strcmp(mode, 'free_fall')
    r_ddot     = generalized_accelerations(1);
    theta_ddot = generalized_accelerations(2);
    psi_ddot   = generalized_accelerations(3);
else
    r_ddot     = 0;   % constrained; r is held fixed by the rolling contact
    theta_ddot = generalized_accelerations(1);
    psi_ddot   = generalized_accelerations(2);
end

dx = [r_dot; r_ddot; theta_dot; theta_ddot; psi_dot; psi_ddot];
end
