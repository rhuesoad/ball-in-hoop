function run = run_scenario(scn, params, cfg)
%RUN_SCENARIO  Shared execution pipeline for every scenario file.
%
%   run = run_scenario(scn)
%   run = run_scenario(scn, params, cfg)
%
%   The only place that knows how scenario declaration -> simulation ->
%   result connects (docs/MODEL.md sec. 9): validates scn against the
%   schema, builds the controller (if any), integrates, assembles the
%   run struct, checks the realised mode sequence against
%   scn.expected_mode_sequence (if declared), saves, and plots.
%
%   Every scenario file is runnable standalone because this function
%   puts the whole project on the MATLAB path itself, from its own
%   location, regardless of the caller's current directory.
%
%   Inputs
%   ------
%   scn    : scenario struct, see docs/MODEL.md sec. 9
%   params : physical parameters (default: ball_hoop_params())
%   cfg    : numerical settings (default: sim_config())
%
%   Outputs
%   -------
%   run : struct obeying results/RUN_SCHEMA.md, extended per
%         docs/MODEL.md sec. 9.7

    project_root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(project_root));

    if nargin < 2 || isempty(params), params = ball_hoop_params(); end
    if nargin < 3 || isempty(cfg),    cfg    = sim_config();       end

    % ball_spin_store.m keeps the ball's spin across a free flight in a
    % persistent variable, which would otherwise survive from one
    % simulation to the next and let a previous run's liftoff set this
    % run's impact.
    ball_spin_store('reset');

    validate_scenario(scn);

    run_cfg = cfg;
    if isfield(scn, 't_end')
        run_cfg.t_end = scn.t_end;
    end

    is_closed_loop = isfield(scn, 'controller') && ~isempty(scn.controller);
    
    if is_closed_loop
        [run, sol] = run_closed_loop(scn, params, run_cfg);
    else
        [run, sol] = run_open_loop(scn, params, run_cfg);
    end

    if isfield(scn, 'expected_mode_sequence')
        run.mode_sequence_check = mode_sequence_check_with_state( ...
            scn.expected_mode_sequence, sol, scn.x0);
        if ~run.mode_sequence_check.matches
            fprintf(['  Mode sequence diverged at interval %d: expected ''%s'', ' ...
                      'got ''%s'' (t = %.4f s).\n'], ...
                run.mode_sequence_check.divergence_index, ...
                run.mode_sequence_check.expected_mode, ...
                run.mode_sequence_check.actual_mode, ...
                run.mode_sequence_check.divergence_time_s);
        end
    end

    save_run(run, fullfile(project_root, 'results', 'runs'));

    if is_closed_loop
        % Tracking/error/torque/phase-portrait summary -- needs a
        % reference, which every closed-loop run has (docs/MODEL.md sec. 9.4).
        plot_run(run, struct('params', params));
    else
        % Open-loop runs have no reference to track (run.ref_psi_rad is
        % all-NaN, results/RUN_SCHEMA.md): the natural view is the ball's
        % path over the hoop geometry instead.
        map2plot = @(x_phys, y_phys) deal(y_phys, -x_phys);
        plot_trajectory(params, run.x(:,1), run.x(:,3), run.x(:,5), sol.Xe, map2plot);
    end
end


function [run, sol] = run_open_loop(scn, params, cfg)
%RUN_OPEN_LOOP  Passive-dynamics integration: no controller, tau = 0.
    tau_fun = @(t, x) 0;
    sol = ball_hoop_ode(scn.x0, tau_fun, @ball_hoop_dynamics, params, cfg, scn.mode);

    n = numel(sol.t);
    tau_log = zeros(n, 1);
    u_log   = zeros(n, 1);

    run = build_run(sol.t, sol.X, tau_log, u_log, scn, params, cfg);
    run.mode_intervals = sol.mode_intervals;
    run.ode_segments   = sol.ode_segments;
end


function [run, sol] = run_closed_loop(scn, params, cfg)
%RUN_CLOSED_LOOP  Joint physical+controller(+motor-lag) integration.
    n_phys = 6;

    [u_fun, ctrl_params, n_ctrl_states, ref, diagnostics] = build_controller(scn, params, cfg);

    % cfg.motor_model = 'first_order' (off by default -- sim_config.m
    % defaults to 'ideal') adds one more integrated state, the lagged
    % motor torque tau_lag (docs/MODEL.md sec. 10.5, motor_inner_loop.m).
    has_motor_lag  = strcmp(cfg.motor_model, 'first_order');
    n_motor_states = double(has_motor_lag);

    % SAMPLED-DATA LOOP. The bench is a discrete controller driving a
    % continuous plant: the camera and encoder deliver one sample every
    % cfg.h_ctrl, the controller emits one torque, and that torque is held
    % while the ball keeps moving. This loop is that, literally -- measure
    % at t_k, compute u_k and tau_k once, then integrate the continuous
    % plant across [t_k, t_k+h] with both frozen.
    %
    % It replaces a zero-order hold that lived inside the ODE right-hand
    % side and kept its tick state in `persistent` variables. That was not
    % merely inelegant: ode15s evaluates the dynamics at trial points and
    % rejects steps, stepping non-monotonically in t, so the persistent
    % tick could be advanced from a trial point that was then thrown away
    % -- leaving the plant driven by a command belonging to no real tick.
    % Freezing the command outside the solver removes the failure mode
    % rather than papering over it.
    %
    % It also makes the logged command the applied command by
    % construction. The previous logging path re-evaluated the controller
    % at every output time, producing a smooth curve where the plant had
    % actually seen a staircase.
    %
    % Mode transitions stay INSIDE a tick: ball_hoop_ode still does its own
    % event detection across each 20 ms interval, because the ball does not
    % wait for the next camera frame to leave the hoop. The controller only
    % learns about it at the following tick -- exactly the bench's
    % situation.
    h      = cfg.h_ctrl;
    x_aug  = [scn.x0(:); zeros(n_ctrl_states, 1); zeros(n_motor_states, 1)];
    mode_k = scn.mode;

    t_all = []; X_aug_all = []; tau_log = []; u_log = [];
    te_all = []; Xe_all = [];
    mode_intervals = struct('t_start', {}, 't_end', {}, 'mode', {});
    ode_segments = {};
    x_history = zeros(0, 1 + n_phys);

    n_ticks = max(1, ceil((cfg.t_end - cfg.t0) / h - 1e-9));

    for k = 0:n_ticks-1
        t_k    = cfg.t0 + k * h;
        t_next = min(t_k + h, cfg.t_end);
        if t_next <= t_k
            break;
        end

        x_phys = x_aug(1:n_phys);
        z      = x_aug(n_phys+1 : n_phys+n_ctrl_states);

        % --- Measure (once per tick, degraded per cfg) and decide ---
        x_history(end+1, :) = [t_k, x_phys.']; %#ok<AGROW>
        x_measured = measured_state_local(x_history, t_k, cfg);
        [u_k, dz_dt_k, ~] = u_fun(t_k, x_measured, z);

        if has_motor_lag
            tau_k = motor_inner_loop(u_k, x_phys, params, cfg, x_aug(end));
        else
            tau_k = motor_inner_loop(u_k, x_phys, params, cfg);
        end

        % --- Hold, and let the continuous plant run for one period ---
        cfg_tick       = cfg;
        cfg_tick.t0    = t_k;
        cfg_tick.t_end = t_next;
        tick_dynamics = @(t, xa, ~, p, m) frozen_dynamics_local( ...
            t, xa, u_k, dz_dt_k, p, cfg, m, n_phys, n_ctrl_states);
        sol_k = ball_hoop_ode(x_aug, @(t, x) tau_k, tick_dynamics, params, cfg_tick, mode_k);

        % Samples of this tick, with the command that actually produced
        % them. sol_k's first sample is t_k itself, already the state we
        % measured, so it is kept here and skipped on later ticks by
        % dropping the previous tick's duplicate endpoint instead.
        n_k = numel(sol_k.t);
        t_all     = [t_all;     sol_k.t(:)];                 %#ok<AGROW>
        X_aug_all = [X_aug_all; sol_k.X];                    %#ok<AGROW>
        u_log     = [u_log;     repmat(u_k,   n_k, 1)];      %#ok<AGROW>
        tau_log   = [tau_log;   repmat(tau_k, n_k, 1)];      %#ok<AGROW>
        te_all    = [te_all;    sol_k.te(:)];                %#ok<AGROW>
        if ~isempty(sol_k.Xe)
            Xe_all = [Xe_all; sol_k.Xe];                     %#ok<AGROW>
        end
        mode_intervals = [mode_intervals, sol_k.mode_intervals];   %#ok<AGROW>
        ode_segments   = [ode_segments,   sol_k.ode_segments];     %#ok<AGROW>

        x_aug  = sol_k.X(end, :).';
        mode_k = sol_k.mode_intervals(end).mode;
    end

    % Drop the duplicated tick boundaries (each tick restarts at the
    % previous tick's final time). Keeping the LAST of each duplicate pair
    % would hide the staircase; keeping the first preserves the held value
    % up to the instant the next command takes effect.
    keep = [true; diff(t_all) > 0];
    t_all     = t_all(keep);
    X_aug_all = X_aug_all(keep, :);
    u_log     = u_log(keep);
    tau_log   = tau_log(keep);

    sol = struct('t', t_all, 'X', X_aug_all, 'te', te_all, 'Xe', Xe_all, ...
        'mode_intervals', merge_mode_intervals_local(mode_intervals), ...
        'ode_segments', {ode_segments});

    X_all = X_aug_all(:, 1:n_phys);
    Z_all = X_aug_all(:, n_phys+1 : n_phys+n_ctrl_states);

    scn_for_run = scn;
    scn_for_run.reference = ref;   % the reference actually used, incl. TVLQR's derived one
    run = build_run(t_all, X_all, tau_log, u_log, scn_for_run, params, cfg);
    run.z_ctrl         = Z_all;
    run.mode_intervals = sol.mode_intervals;
    % NOTE: unlike run.x, deval(run.ode_segments{k}, tq) returns the full
    % AUGMENTED state [x_phys; z; tau_lag(if any)] (ball_hoop_ode.m was
    % called with x0_aug) -- take only the first n_phys=6 rows for the
    % physical state.
    run.ode_segments   = sol.ode_segments;

    if has_motor_lag
        run.tau_lag_Nm = X_aug_all(:, n_phys+n_ctrl_states+1);
    end

    % The synthesised controller itself, not just its closed-loop poles.
    % Without this the gains a run was produced with are unrecoverable
    % after the fact: build_controller computes them, run_scenario used to
    % discard them, and the scenario files do not expose scn. They are
    % provenance -- and they are what has to be transcribed to the
    % hardware controller (analysis/verification/extract_gains.m).
    run.ctrl_params = ctrl_params;

    if isfield(diagnostics, 'eig_cl')
        run.eig_cl = diagnostics.eig_cl;
    end
    if isfield(diagnostics, 'solver_info')
        run.solver_info = diagnostics.solver_info;
    end
end


function result = mode_sequence_check_with_state(expected_sequence, sol, x0)
%MODE_SEQUENCE_CHECK_WITH_STATE  Wraps check_mode_sequence.m, filling in
%                                 the state at the divergence point from
%                                 sol (which check_mode_sequence.m itself
%                                 does not see -- docs/MODEL.md sec. 9.5).
    result = check_mode_sequence(expected_sequence, sol.mode_intervals);
    if result.matches
        return;
    end

    if result.divergence_index == 1
        % First interval: the run started already off the expected script.
        result.divergence_state = x0;
    else
        % t_start of interval k (k>1) is, by construction, the (k-1)-th
        % transition time in sol.te -- ball_hoop_ode.m closes an interval
        % and opens the next at exactly that value.
        event_idx = result.divergence_index - 1;
        if event_idx <= size(sol.Xe, 1)
            result.divergence_state = sol.Xe(event_idx, 1:6).';
        else
            % The run produced FEWER transitions than expected -- it ran
            % out of time (or stalled) before reaching the interval that
            % diverged, so no event state exists to report. Diverging off
            % the end of the sequence is a legitimate outcome to describe,
            % not an indexing error.
            result.divergence_state = [];
        end
    end
end


function dx_aug = augmented_dynamics_local(t, x_aug, u_fun, params, cfg, current_mode, n_phys, n_ctrl_states)
%AUGMENTED_DYNAMICS_LOCAL  Joint derivative of (physical state, controller
%                          state, motor-lag state if any), zero-order-hold
%                          outer loop. Adapted from ball_hoop_closed_loop.m's
%                          local function of the same purpose, for the new
%                          scn schema's controller-building path
%                          (build_controller.m).
%
%   Inputs
%   ------
%   t            : current time [s], or the string 'reset' to clear the cache
%   x_aug        : augmented state [x_phys; z; tau_lag(if any)]
%                  (n_phys+n_ctrl_states+n_motor_states x 1)
%   u_fun        : controller handle [u, dz_dt, info] = u_fun(t, x_phys, z)
%   params       : physical parameters
%   cfg          : numerical settings, incl. cfg.h_ctrl (sample period [s])
%                  and cfg.motor_model ('ideal' | 'first_order')
%   current_mode : 'rolling_out' | 'rolling_in_outside' | 'rolling_in_inside' | 'free_fall'
%   n_phys       : number of physical states (6)
%   n_ctrl_states : number of controller internal states
%
%   Outputs
%   -------
%   dx_aug : time derivative of the augmented state, same layout as x_aug

    persistent last_t_k last_u last_dz_dt x_history

    if ischar(t) && strcmp(t, 'reset')
        last_t_k = [];
        last_u = [];
        last_dz_dt = [];
        x_history = [];
        dx_aug = [];
        return;
    end

    x_phys = x_aug(1:n_phys);
    z      = x_aug(n_phys+1 : n_phys+n_ctrl_states);

    t_k = floor(t / cfg.h_ctrl) * cfg.h_ctrl;

    if isempty(last_t_k) || t_k ~= last_t_k
        % Sample the measurement once per new tick (docs/VALIDATION.md) --
        % same reason u_fun itself is only re-evaluated once per tick,
        % see below: ode15s calls this function many times per accepted
        % step during adaptive stepping, and a real sensor doesn't.
        x_history(end+1, :) = [t_k, x_phys.']; %#ok<AGROW>
        x_measured = measured_state_local(x_history, t_k, cfg);
        [last_u, last_dz_dt, ~] = u_fun(t, x_measured, z);
        last_t_k = t_k;
    end
    u     = last_u;
    dz_dt = last_dz_dt;

    if strcmp(cfg.motor_model, 'first_order')
        tau_lag = x_aug(n_phys+n_ctrl_states+1);
        [tau, dtau_lag_dt] = motor_inner_loop(u, x_phys, params, cfg, tau_lag);
        motor_state_dot = dtau_lag_dt;
    else
        tau = motor_inner_loop(u, x_phys, params, cfg);
        motor_state_dot = [];
    end

    dx_phys = ball_hoop_dynamics(t, x_phys, @(tt, xx) tau, params, current_mode);
    dx_aug  = [dx_phys; dz_dt; motor_state_dot];
end


function dx_aug = frozen_dynamics_local(t, x_aug, u_k, dz_dt_k, params, cfg, current_mode, n_phys, n_ctrl_states)
%FROZEN_DYNAMICS_LOCAL  Augmented dynamics over one control period, with
%                       the controller output held constant.
%
%   The counterpart of the sampled-data loop in run_closed_loop: u_k and
%   dz_dt_k were computed once, at the tick, from the measured state, and
%   do not change until the next tick. Nothing here calls the controller,
%   which is the whole point -- the solver may evaluate this function
%   anywhere in the interval, including at trial points it later rejects,
%   without any of that reaching the control law.
    x_phys = x_aug(1:n_phys);

    if strcmp(cfg.motor_model, 'first_order')
        tau_lag = x_aug(n_phys+n_ctrl_states+1);
        [tau, dtau_lag_dt] = motor_inner_loop(u_k, x_phys, params, cfg, tau_lag);
        motor_state_dot = dtau_lag_dt;
    else
        tau = motor_inner_loop(u_k, x_phys, params, cfg);
        motor_state_dot = [];
    end

    dx_phys = ball_hoop_dynamics(t, x_phys, @(tt, xx) tau, params, current_mode);
    dx_aug  = [dx_phys; dz_dt_k; motor_state_dot];
end


function merged = merge_mode_intervals_local(intervals)
%MERGE_MODE_INTERVALS_LOCAL  Joins intervals split only by a tick boundary.
%
%   Each control period is integrated by its own ball_hoop_ode call, so a
%   stretch of unchanging mode arrives as one interval per tick. Those are
%   an artefact of how the loop is driven, not mode history: merge runs of
%   the same mode back into single intervals so mode_intervals still means
%   what the rest of the codebase expects.
    merged = struct('t_start', {}, 't_end', {}, 'mode', {});
    for i = 1:numel(intervals)
        if ~isempty(merged) && strcmp(merged(end).mode, intervals(i).mode)
            merged(end).t_end = intervals(i).t_end;
        else
            merged(end+1) = intervals(i); %#ok<AGROW>
        end
    end
end


function x_measured = measured_state_local(x_history, t_k, cfg)
%MEASURED_STATE_LOCAL  Degrades the true, per-tick-sampled physical state
%                      into what the controller actually sees: delay,
%                      then quantization, then noise (docs/VALIDATION.md).
%                      With every cfg field at its default this returns
%                      x_phys completely unchanged.
%
%   Inputs
%   ------
%   x_history : [t_k, x_phys.'] rows accumulated so far, one per tick,
%               in increasing t_k order (by construction: the caller
%               only appends at a new tick)
%   t_k       : current tick time [s]
%   cfg       : numerical settings, incl. .measurement_delay [s],
%               .theta_quantization_rad, .psi_quantization_rad [rad],
%               .measurement_noise_std (struct, [rad] and [rad/s])
%
%   Outputs
%   -------
%   x_measured : degraded state (n_phys x 1)

    % --- Delay: use the sample from (at least) measurement_delay ago ---
    % Discrete, tick-resolution delay -- consistent with the fact that
    % the controller only ever sees one sample per tick anyway, so a
    % continuous-time interpolated delay would imply sensor resolution
    % the rest of this model doesn't have.
    if cfg.measurement_delay > 0
        target_t = t_k - cfg.measurement_delay;
        idx = find(x_history(:,1) <= target_t + 1e-9, 1, 'last');
        if isempty(idx), idx = 1; end
        x_measured = x_history(idx, 2:end).';
    else
        x_measured = x_history(end, 2:end).';
    end

    % --- Quantization ---
    if cfg.theta_quantization_rad > 0
        x_measured(3) = round(x_measured(3) / cfg.theta_quantization_rad) * cfg.theta_quantization_rad;
    end
    if cfg.psi_quantization_rad > 0
        x_measured(5) = round(x_measured(5) / cfg.psi_quantization_rad) * cfg.psi_quantization_rad;
    end

    % --- Noise (zero-mean Gaussian; see sim_config.m for reproducibility note) ---
    ns = cfg.measurement_noise_std;
    if ns.theta_rad > 0,       x_measured(3) = x_measured(3) + ns.theta_rad       * randn(); end
    if ns.theta_dot_rad_s > 0, x_measured(4) = x_measured(4) + ns.theta_dot_rad_s * randn(); end
    if ns.psi_rad > 0,         x_measured(5) = x_measured(5) + ns.psi_rad         * randn(); end
    if ns.psi_dot_rad_s > 0,   x_measured(6) = x_measured(6) + ns.psi_dot_rad_s   * randn(); end
end


function [tau, u] = call_controller_tau_local(t, x_aug, u_fun, n_phys, n_ctrl_states, params, cfg)
%CALL_CONTROLLER_TAU_LOCAL  Recomputes tau (and the commanded u) at an
%                           accepted solution point, for logging only.
%                           See augmented_dynamics_local's docstring for
%                           why this is a re-evaluation, not a replay of
%                           the zero-order-hold cache. When any
%                           measurement-realism field is active
%                           (sim_config.m), this is an approximation on
%                           top of that existing one: quantization and
%                           noise are reapplied fresh (noise draws a new
%                           sample, so the logged tau/u track the noise
%                           statistics but not the exact realised
%                           sequence), and measurement_delay is NOT
%                           reapplied at all (this function only sees one
%                           instant, not the history a delay needs) --
%                           the logged tau/u under-report the effect of a
%                           nonzero delay. Exact reproduction of what was
%                           actually applied during integration would
%                           require logging it live, not recomputing it.
%
%   Inputs
%   ------
%   t      : current time [s]
%   x_aug  : augmented state [x_phys; z; tau_lag(if any)]
%   u_fun  : controller handle [u, dz_dt, info] = u_fun(t, x_phys, z)
%   n_phys : number of physical states (6)
%   n_ctrl_states : number of controller internal states
%   params : physical parameters
%   cfg    : numerical settings, incl. cfg.tau_max [N.m], cfg.motor_model
%
%   Outputs
%   -------
%   tau : applied (saturated) motor torque [N.m]
%   u   : commanded angular acceleration [rad/s^2]
    x_phys = x_aug(1:n_phys);
    z      = x_aug(n_phys+1 : n_phys+n_ctrl_states);

    x_measured = measured_state_local([t, x_phys.'], t, cfg);
    [u, ~, ~]  = u_fun(t, x_measured, z);

    if strcmp(cfg.motor_model, 'first_order')
        tau_lag = x_aug(n_phys+n_ctrl_states+1);
        tau     = motor_inner_loop(u, x_phys, params, cfg, tau_lag);
    else
        tau = motor_inner_loop(u, x_phys, params, cfg);
    end
end
