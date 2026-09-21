function [u, dz_dt, info] = t4_flying_ball_controller(t, x_phys, ctrl)
    r     = x_phys(1);
    r_dot = x_phys(2);
    tol   = ctrl.contact_tol;

    in_contact = abs(r_dot) < tol;

    if in_contact && abs(r - ctrl.R_out) < tol && t <= ctrl.ref_traj.Tf
        % --- Phase 1: rolling on the outer hoop, tracking the planned
        %     energy-injection trajectory. It deliberately STOPS short of
        %     liftoff: past Tf, the ball is sill on the hoop but coasting,
        %     so the controller would hold the terminal gain and fight the 
        %     ball for overshooting. But in this case, it is MEANT to do it
        [u, dz_dt, info] = tvlqr_controller(t, x_phys, ctrl.ref_traj, ctrl.ctrl_params);
        info.phase = 1;

    elseif in_contact && abs(r - ctrl.R_in_i) < tol
        % --- Phase 3: caught. Stabilise at the bottom of the inner
        %     hoop inside surface, task = T1 variant. ---
        cf_active = isfield(ctrl, 'catch_follow') && ~isempty(ctrl.catch_follow) ...
                    && abs(x_phys(6)) > ctrl.catch_follow.psi_dot_switch;
        if cf_active
            % --- Phase 3a: RUN WITH THE BALL, SLIGHTLY SLOWER. ---
            % The hoop first follows the ball to avoid that braking breaks
            % the trajectory due to the ball slipping because it has too
            % much speed. The almost-matching motion allows the ball to
            % dissipate some energy through rolling contact. Once the ball
            % is slow enough, phase 3b takes place (actual regulation).
            cf = ctrl.catch_follow;
            theta_dot     = x_phys(4);
            theta_dot_ref = cf.ratio * x_phys(6);

            u_hi = ( ctrl.sat.tau_max - ctrl.sat.motor_friction * theta_dot) / ctrl.sat.total_inertia;
            u_lo = (-ctrl.sat.tau_max - ctrl.sat.motor_friction * theta_dot) / ctrl.sat.total_inertia;
            u_unsat = -cf.kd * (theta_dot - theta_dot_ref);
            u = max(u_lo, min(u_hi, u_unsat));

            dz_dt = zeros(0, 1);
            info = struct('phase', 3, 'sub', 'follow', 'u', u, ...
                          'u_unsat', u_unsat, 'saturated', u ~= u_unsat, ...
                          'theta_dot_ref', theta_dot_ref);
        else
            % --- Phase 3b: regulate at the bottom. ---
            % Regulator just like T1. 
            theta_park = ctrl.theta_park_offset ...
                       + 2*pi * round((x_phys(3) - ctrl.theta_park_offset) / (2*pi));

            x_parked    = x_phys;
            x_parked(3) = x_phys(3) - theta_park;   % error folded into (-pi, pi]

            [u, dz_dt, info] = lqr_controller(t, x_parked, ctrl.ref_zero, ctrl.lqr);
            info.phase      = 3;
            info.sub        = 'regulate';
            info.theta_park = theta_park;
        end

    else
        % --- Phase 2: free flight. ---
        theta     = x_phys(3);
        theta_dot = x_phys(4);

        % Saturation
        u_hi = ( ctrl.sat.tau_max - ctrl.sat.motor_friction * theta_dot) / ctrl.sat.total_inertia;
        u_lo = (-ctrl.sat.tau_max - ctrl.sat.motor_friction * theta_dot) / ctrl.sat.total_inertia;
        u_bang = min(abs(u_hi), abs(u_lo));

        if isfield(ctrl, 'theta_dot_hold') && ~isempty(ctrl.theta_dot_hold)
            % Speed holds. Bringing the hoop to rest during flight is the
            % wrong target to reach: the impact gives back energy at
            % impact, and the only thing that cancels it if the ball must
            % brake, is the rotating hoop itself. Stopping the hoops during
            % the flight makes it unable to dissipate enough energy to
            % cancel the ball's motion. It would therefore be caught and
            % thrown back out through the gap. 
            % 
            % The real problem of free-flight is therefore holding the
            % speed. Phase 1 already picked the release angle inside a
            % window that makes landing admissible. 
            t_flight = t - ctrl.ref_traj.Tf;
            theta_ref = ctrl.theta_hold + ctrl.theta_dot_hold * t_flight;
            theta_error = wrapToPi(theta - theta_ref);
            u_unsat = -ctrl.kp_hold * theta_error ...
                      - ctrl.kd_hold * (theta_dot - ctrl.theta_dot_hold);
            s = NaN;
        else
            theta_error = wrapToPi(theta - ctrl.theta_hold);
            s = theta_error + theta_dot * abs(theta_dot) / (2 * u_bang);

            if abs(s) > ctrl.slew_boundary
                u_unsat = -u_bang * sign(s);
            else
                u_unsat = -ctrl.kp_hold * theta_error - ctrl.kd_hold * theta_dot;
            end
        end
        u = max(u_lo, min(u_hi, u_unsat));

        dz_dt = zeros(0, 1);
        info = struct('phase', 2, 'u_unsat', u_unsat, 'u', u, ...
                      'saturated', u ~= u_unsat, 'theta_error', theta_error, ...
                      'switching_function', s);
    end
end
