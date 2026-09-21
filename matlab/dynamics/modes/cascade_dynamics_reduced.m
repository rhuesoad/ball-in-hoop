function dx = cascade_dynamics_reduced(x, u, params, mode)
%CASCADE_DYNAMICS_REDUCED  Reduced (theta, theta_dot, psi, psi_dot)
%                          cascade dynamics used by trajectory collocation.
%
%   dx = cascade_dynamics_reduced(x, u, params, mode)
%
%   Inputs
%   ------
%   x      : reduced state [theta; theta_dot; psi; psi_dot]
%   u      : control input, u = theta_ddot                      [rad/s^2]
%   params : physical parameters struct (see ball_hoop_params)
%   mode   : (optional) 'rolling_out' | 'rolling_in_outside' | 'rolling_in_inside'
%            Default 'rolling_out', preserving every call site that
%            existed before this parameter was added (all of which were
%            hardcoded to the outer surface -- see docs/AUDIT.md item C7).
%
%   Output
%   ------
%   dx : [theta_dot; u; psi_dot; psi_ddot]

    if nargin < 4
        mode = 'rolling_out';
    end

    theta_dot = x(2);
    psi       = x(3);
    psi_dot   = x(4);

    [M, C, G] = rolling_matrices(params, psi, mode);

    psi_ddot = -(M(2,1)*u + C(2,1)*theta_dot + C(2,2)*psi_dot + G(2)) / M(2,2);
    dx = [theta_dot; u; psi_dot; psi_ddot];

end
