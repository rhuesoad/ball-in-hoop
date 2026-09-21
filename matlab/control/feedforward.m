function u_ff = feedforward(psi_ref, psi_dot_ref, psi_ddot_ref, ...
                                  theta_dot, ctrl_params)
%FEEDFORWARD_EXACT  Open-loop acceleration that realizes a psi trajectory.
%
%   u_ff = feedforward_exact(psi_ref, psi_dot_ref, psi_ddot_ref, ...
%                            theta_dot, ctrl_params)
%
%   Exact inversion of the ball equation in cascade form. The reduced
%   (theta, psi) dynamics on a rolling surface read:
%
%       M(2,1)*theta_ddot + M(2,2)*psi_ddot + C(2,1)*theta_dot
%                         + C(2,2)*psi_dot + m*g*R_eff*sin(psi) = 0
%
%   Imposing the desired trajectory (psi, psi_dot, psi_ddot) =
%   (psi_ref, psi_dot_ref, psi_ddot_ref) and the input u = theta_ddot,
%   we solve for u:
%
%       u_ff = -[ M(2,2)*psi_ddot_ref + C(2,1)*theta_dot
%                 + C(2,2)*psi_dot_ref + m*g*R_eff*sin(psi_ref) ] / M(2,1)
%
%   The Coriolis term C(2,1)*theta_dot uses the LIVE hoop velocity
%   theta_dot. There is no closed form for theta_dot along the trajectory
%   (it is itself driven by u_ff), so it is read from the current state at
%   each call. Evaluated on the true trajectory, u_ff reproduces psi_ref in
%   open loop up to the actuator saturation; the LQR corrects the residual.
%
%   Inputs
%   ------
%   psi_ref      : desired ball angle      [rad]
%   psi_dot_ref  : desired ball angular rate [rad/s]
%   psi_ddot_ref : desired ball angular accel [rad/s^2]
%   theta_dot    : live hoop angular velocity [rad/s]
%   ctrl_params  : must contain .mode, .R_eff [m], .params
%
%   Output
%   ------
%   u_ff : feedforward angular acceleration u = theta_ddot_ref [rad/s^2]
%
%   Reference: standard computed-torque / inverse-dynamics feedforward,
%   e.g. Lynch & Park (2017), Modern Robotics, sec. 8.4 (feedforward
%   linearizing control) and 11.4 (feedforward plus feedback).

    params = ctrl_params.params;
    mode   = ctrl_params.mode;
    R_eff  = ctrl_params.R_eff;

    [M, C, ~] = rolling_matrices(params, psi_ref, mode);

    u_ff = -( M(2,2) * psi_ddot_ref ...
            + C(2,1) * theta_dot ...
            + C(2,2) * psi_dot_ref ...
            + params.ball_mass * params.gravity * R_eff * sin(psi_ref) ) ...
            / M(2,1);
end