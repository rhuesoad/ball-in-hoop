function cfg = sim_config()
%SIM_CONFIG  Numerical settings for the ball-hoop simulation.
%
%   cfg = sim_config()
%
%   Inputs:  none
%   Outputs: cfg (struct) — all numerical settings in SI units
%
%   Contains everything that is NOT physical: time horizon, actuator limits,
%   ODE solver tolerances, and event-detection constants used by the hybrid
%   integrator. No physical quantities belong here.
%
%   Separation rationale: the same physical model (ball_hoop_params) can be
%   re-simulated under a different numerical regime — finer grid, tighter
%   tolerances, longer horizon — by editing only this file.

    %% --- Simulation time horizon ---
    cfg.t0    = 0;                                          % [s]  start time
    cfg.t_end = 10.0;                                        % [s]  end time

    %% --- Actuator saturation ---
    % Physical limit of the motor + ODrive driver combination. Every
    % controller reads this value from cfg; never hard-code it elsewhere.
    %
    % tau_max = Kt * I_configured, with:
    %   Kt = 0.031 N.m/A -- ball_hoop_params.m's motor_torque_constant
    %        (D5065 270KV). Kept as a literal here rather than read from
    %        params: sim_config.m is deliberately independent of
    %        ball_hoop_params.m (see the separation rationale above). Keep
    %        the two in sync by hand if Kt is ever re-identified.
    %   I  = 17.5 A -- planned ODrive S1 current limit once both 24V/10A
    %        supplies are wired in. PROVISIONAL: the driver is currently
    %        configured at 8A on a single supply, which would give
    %        tau_max = 0.248 N.m instead. Update this value once the
    %        supplies are actually wired and the current limit raised, and
    %        confirm the two supplies really share current before relying
    %        on 17.5 A.
    cfg.tau_max = 0.310;                                   % [N.m]
    
    %% --- Command and speed limits ---
    cfg.theta_dot_max = 35.0;    % [rad/s] maximum hoop speed
    cfg.u_slip_max = 400.0;           % [rad/s^2] no-slip / command acceleration limit
    
    %% --- Motor inner loop ---
    % Cascade architecture: outer controller produces theta_ddot_ref,
    % converted to torque by an emulated inner velocity loop.
    cfg.motor_model = 'first_order';        % 'ideal' | 'first_order'
    cfg.T_loop      = 0.015;           % [s]  exp., inner-loop time constant
                                       %     (used only by 'first_order' model;
                                       %      was 0.010, an assumed value)
    
    %% --- ODE solver ---
    % ode15s is required: the rolling-contact DAE structure makes the system
    % stiff when the ball transitions between modes.
    cfg.solver = 'ode15s';
    cfg.RelTol = 1e-5;
    cfg.AbsTol = 1e-6;

    %% --- Hybrid event-detection constants ---

    % After a mode transition, restart the integrator at t_event + SMALL_DT
    % so the same zero-crossing is not immediately re-detected in the new mode.
    cfg.SMALL_DT = 1e-8;                                    % [s]

    % Small displacement added to the state at a transition so the new mode's
    % contact condition is strictly satisfied from the first step, preventing
    % spurious immediate re-entry into the previous mode.
    cfg.TOL_NUDGE = 1e-8;                                   % [m]

    % If an event is detected at time t0 itself (degenerate crossing), the
    % integrator skips it. This threshold distinguishes a real crossing from
    % floating-point noise at the segment boundary.
    cfg.TOL_DEGENERATE = 1e-12;                             % [s]

    % Safety cap on the total number of integration segments. Each mode switch
    % opens a new segment, so a bug that cycles between two modes would
    % otherwise run forever.
    % Note for Python port: treat as integer (use int, not float).
    cfg.MAX_ITER = 1e6;                                     % [dimensionless]

    %% --- Impact and event-detection tolerances ---

    % Angular padding added to both edges of the hole boundary in detect_events
    % so a ball approaching the edge is not caught exactly at the geometric boundary.
    cfg.ANGLE_PAD = 1e-3;                                   % [rad]

    % Tolerance for the hole-boundary test inside handle_impact; guards against
    % floating-point noise when psi is numerically exactly at hole_start/end.
    cfg.TOL_ANGLE = 1e-6;                                   % [rad]

    % Tolerance for the contact-radius test in handle_impact; the ball's
    % radial position at an event is not exactly on the hoop surface due to
    % solver step size, so a small band is needed.
    cfg.TOL_RADIUS = 1e-4;                                  % [m]

    %% --- Sampling ---
    cfg.h_ctrl = 1/50;   % [s]  outer control loop sample period (50 Hz)

    %% --- Measurement realism (Phase 6 -- docs/VALIDATION.md) ---
    % All off/neutral by default: with every field at its default below,
    % the controller sees the true state exactly, so the closed loop is
    % unaffected by anything in this block. Set any of them for a
    % "realistic simulation" comparison run (docs/VALIDATION.md);
    % run_scenario.m's closed-loop path is the only place that reads
    % them, applying each to the state *fed to the controller*, never to
    % the physical dynamics themselves -- the ball still falls under its
    % real, uncorrupted state, only the controller's sensing is degraded.

    % Camera/encoder latency: the controller acts on the state from
    % measurement_delay seconds ago rather than the current one. 0 = no
    % delay.
    %
    % 19.3 ms is the measured camera pipeline latency, replacing the round
    % 20 ms that stood here. RESOLUTION: measured_state_local (run_scenario.m)
    % resolves the delay to whole ticks, so any value in (0, h_ctrl] behaves
    % as exactly one tick -- 19.3 ms and 20 ms produce bit-identical runs at
    % the default h_ctrl = 20 ms. The value is set honestly here so the
    % provenance is right, not because the sim can currently tell them apart.
    cfg.measurement_delay = 0.0193;                             % [s]  exp.

    % Encoder/vision angular resolution: the measured angle is rounded to
    % the nearest multiple of this step before reaching the controller.
    % 0 = unquantized (infinite resolution). Values must come from the
    % actual hardware (encoder counts/rev, camera pixels/rev at the
    % working distance) -- do not guess them.
    cfg.theta_quantization_rad = 0;                         % [rad]
    cfg.psi_quantization_rad   = 0;                         % [rad]

    % Measurement noise (zero-mean Gaussian, added to the state the
    % controller sees). All zero = no noise. Variances must come from
    % calib_analyze.py's actual characterisation of the sensors, not an
    % assumed number (docs/VALIDATION.md) -- MATLAB's rng state is not
    % managed here, so callers wanting reproducible noise must rng(seed)
    % themselves before calling run_scenario.
    cfg.measurement_noise_std.theta_rad       = 0;          % [rad]
    cfg.measurement_noise_std.theta_dot_rad_s = 0;          % [rad/s]
    cfg.measurement_noise_std.psi_rad         = 0;          % [rad]
    cfg.measurement_noise_std.psi_dot_rad_s   = 0;          % [rad/s]
end
