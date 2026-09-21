function [tau, dtau_lag_dt] = motor_inner_loop(u, x_phys, params, cfg, tau_lag)
%MOTOR_INNER_LOOP  Converts an angular acceleration setpoint into a motor
%                  torque, emulating the inner velocity loop of the ODrive.
%
%   tau = motor_inner_loop(u, x_phys, params, cfg)
%   [tau, dtau_lag_dt] = motor_inner_loop(u, x_phys, params, cfg, tau_lag)
%
%   Inputs
%   ------
%   u       :  angular acceleration setpoint    [rad/s^2]
%             (output of the outer controller PID/LQR/TVLQR)
%   x_phys  : full physical state
%             [r; r_dot; theta; theta_dot; psi; psi_dot]
%   params  : physical parameters from ball_hoop_params()
%   cfg     : numerical settings struct, must contain:
%               cfg.motor_model : 'ideal' | 'first_order'
%               cfg.T_loop      : inner-loop time constant [s]
%                                 (only used by 'first_order')
%   tau_lag : current value of the lagged-torque state [N.m] --
%             required for 'first_order', ignored for 'ideal'. This is
%             an extra state the caller must integrate (see below); it
%             is NOT computed internally, because this function is
%             called fresh at every RHS evaluation with no memory of
%             its own between calls.
%
%   Outputs
%   -------
%   tau         : motor torque to apply to the dynamics, saturated to
%                 +/- cfg.tau_max                                 [N.m]
%   dtau_lag_dt : time derivative of the tau_lag state             [N.m/s]
%                 -- 'first_order' only; callers using 'ideal' have no
%                 extra state to integrate and can ignore this output.
%
%   This is the only place in the codebase that converts u -> tau, and
%   the only place that saturates it (Phase 2 architecture target: the
%   closed-loop pipeline used to have two call sites that each
%   saturated separately after calling this function -- identical
%   logic, duplicated).
%
%   Theory
%   ------
%   In the cascade architecture (Gurtner & Zemanek 2017), the outer
%   controller produces an acceleration setpoint u = theta_ddot_ref.
%   This setpoint must be converted into a torque before being passed to
%   the full Lagrangian dynamics (which expects tau as input).
%
%   On the real hardware, this conversion is performed by the ODrive's
%   internal velocity loop. In simulation, we emulate this conversion
%   here. Two models are available:
%
%   - 'ideal' :
%       Assumes the inner loop is infinitely fast (perfect timescale
%       separation). The torque is computed by algebraic inversion of
%       the simplified motor equation:
%           I_total * theta_ddot + b_motor * theta_dot = tau
%       which gives
%           tau = I_total * u + b_motor * theta_dot
%       This ignores the coupling terms with the ball, just as the
%       ODrive's velocity loop would (it does not know about the ball).
%
%   - 'first_order' :
%       Models the inner loop as a first-order lag with time constant
%       T_loop, capturing imperfect timescale separation. The torque
%       the loop is instantaneously trying to reach is the same
%       tau_ideal the 'ideal' model outputs; the torque actually
%       delivered, tau_lag, tracks it through
%           T_loop * d(tau_lag)/dt + tau_lag = tau_ideal(t)
%       so tau_lag is a state the caller must add to the integrated
%       system (docs/MODEL.md sec. 10.5) -- this function only supplies
%       tau_ideal and dtau_lag_dt = (tau_ideal - tau_lag)/T_loop; it
%       does not integrate anything itself. The delivered torque
%       (saturated, applied to the physical dynamics) is tau_lag
%       itself, not tau_ideal. Used in a validation study to quantify
%       the impact of finite loop bandwidth -- off by default
%       (cfg.motor_model = 'ideal').
%
%   Reference: Khalil (2002), Nonlinear Systems, ch. 11, on timescale
%   separation and singular perturbation methods.

    theta_dot = x_phys(4);

    % Algebraic inversion of the simplified motor equation (both models
    % use this as the loop's instantaneous target torque).
    tau_ideal = params.total_inertia * u + params.motor_friction * theta_dot;

    switch cfg.motor_model
        case 'ideal'
            tau = tau_ideal;
            dtau_lag_dt = [];

        case 'first_order'
            if nargin < 5 || isempty(tau_lag)
                error(['motor_inner_loop: cfg.motor_model = ''first_order'' requires ' ...
                       'the tau_lag state as a 5th argument -- the caller must integrate ' ...
                       'it as an extra state (docs/MODEL.md sec. 10.5).']);
            end
            dtau_lag_dt = (tau_ideal - tau_lag) / cfg.T_loop;
            tau = tau_lag;

        otherwise
            error(['motor_inner_loop: unknown motor_model ''%s''. ' ...
                  'Valid options are ''ideal'' and ''first_order''.'], ...
                  cfg.motor_model);
    end

    % Actuator saturation -- the only place this happens (see docstring).
    % Applied to whichever torque is about to drive the dynamics: tau_ideal
    % itself for 'ideal', the lagged (already-delayed) torque for
    % 'first_order' -- never to tau_ideal in the 'first_order' case, which
    % would saturate a value that's never actually delivered.
    tau = max(-cfg.tau_max, min(cfg.tau_max, tau));
end