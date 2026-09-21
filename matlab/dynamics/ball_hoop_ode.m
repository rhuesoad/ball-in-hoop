function sol = ball_hoop_ode( ...
        x0, tau_fun, dynamics_fun, params, cfg, initial_mode)
%BALL_HOOP_ODE  Integrates the ball-hoop system with hybrid mode transitions.
%
%   Inputs
%   ------
%   x0            : initial state vector
%   tau_fun       : function handle  tau_fun(t, x) -> torque [N.m]
%                   In closed loop, x may be an augmented state (physical +
%                   controller). In open loop, x is the 6-state physical
%                   vector and tau_fun typically ignores x.
%   dynamics_fun  : function handle  dx_dt = dynamics_fun(t, x, tau_fun, params, mode)
%                   Returns the derivative of x. Lets the caller decide
%                   whether x is purely physical (open loop) or augmented
%                   with controller states (closed loop).
%   params        : physical parameters struct
%   cfg           : numerical settings struct
%   initial_mode  : starting dynamic mode — one of:
%                     'free_fall' | 'rolling_out' | 'rolling_in_outside' | 'rolling_in_inside'
%
%   Output
%   ------
%   sol : struct with fields
%           .t              : time vector [s]                    (N x 1)
%           .X              : state matrix                       (N x size(x0,1))
%           .te             : times at which mode transitions occurred [s]
%           .Xe             : states at mode transitions          (M x size(x0,1))
%           .mode_intervals : mode history as a struct array, one entry
%                              per contiguous run of a single mode:
%                                .t_start, .t_end : interval bounds [s],
%                                                    set at the exact
%                                                    transition times, so
%                                                    consecutive intervals
%                                                    tile [t(1), t(end)]
%                                                    with no gaps/overlaps
%                                .mode            : mode string active
%                                                    throughout the interval
%                              (Phase 2 architecture target; previously
%                              a per-sample cell array, docs/AUDIT.md
%                              item A4)
%           .ode_segments   : cell array of continuous ODE solution
%                              structs, one per raw solver call (i.e. one
%                              per while-loop iteration below; a run with
%                              no mode transitions has exactly one).
%                              Each is deval-compatible: deval(sol,tq)
%                              interpolates with the solver's own dense
%                              output, not by fitting a curve through the
%                              discrete .t/.X samples. Doesn't
%                              necessarily line up 1:1 with
%                              .mode_intervals (a degenerate-event retry,
%                              below, stays in the same mode but still
%                              starts a new segment); consumers that need
%                              "the segment covering time t" should
%                              search by each segment's own [x(1) x(end)]
%                              range, not by index (Phase 6 -- see
%                              docs/MODEL.md sec. 10.4).

    t_all  = [];
    X_all  = [];
    te_all = [];
    Xe_all = [];
    ode_segments = {};

    mode_intervals = struct('t_start', {}, 't_end', {}, 'mode', {});
    interval_start = cfg.t0;
    interval_mode  = initial_mode;

    t_current    = cfg.t0;
    t_end        = cfg.t_end;
    current_mode = initial_mode;
    x_current    = x0;
    iter         = 0;

    base_opts = odeset('RelTol', cfg.RelTol, 'AbsTol', cfg.AbsTol);

    while t_current < t_end && iter < cfg.MAX_ITER
        iter = iter + 1;
        previous_mode = current_mode;

        seg_opts = odeset(base_opts, ...
            'Events', @(t, x) detect_events(t, x, params, cfg, current_mode));

        % Single-output ("struct") form: same solve as the previous
        % [t,X,te,Xe,ie] call, but the returned struct is also usable
        % directly with deval (sol.ode_segments, Phase 6) -- MATLAB's ODE
        % solvers only attach the dense-output data needed for deval when
        % called this way, not when called with explicit output arguments.
        sol_seg = feval(cfg.solver, ...
            @(t, x) dynamics_fun(t, x, tau_fun, params, current_mode), ...
            [t_current, t_end], x_current, seg_opts);
        ode_segments{end+1} = sol_seg; %#ok<AGROW>

        t_seg = sol_seg.x.';
        X_seg = sol_seg.y.';
        te    = sol_seg.xe.';
        Xe    = sol_seg.ye.';

        if ~isempty(t_all) && numel(t_seg) > 1
            t_seg = t_seg(2:end);
            X_seg = X_seg(2:end, :);
        end
        if ~isempty(t_seg)
            t_all = [t_all; t_seg];
            X_all = [X_all; X_seg];
        end

        if isempty(te)
            break;
        end

        if abs(te(end) - t_current) < cfg.TOL_DEGENERATE
            t_current = te(end) + cfg.SMALL_DT;
            x_current = Xe(end, :)';
            continue;
        end

        [current_mode, x_current] = transition_state( ...
            previous_mode, Xe(end, :)', params, cfg);

        if ~strcmp(previous_mode, current_mode)
            te_all = [te_all; te];
            Xe_all = [Xe_all; Xe];

            % Close the interval that just ended, exactly at the
            % transition time, so the next interval starts at the same
            % value -- guarantees no gap/overlap between intervals.
            mode_intervals(end+1) = struct( ...
                't_start', interval_start, 't_end', te(end), 'mode', interval_mode); %#ok<AGROW>
            interval_start = te(end);
            interval_mode  = current_mode;
        end

        t_current = te(end) + cfg.SMALL_DT;

    end

    if iter >= cfg.MAX_ITER
        warning(['ball_hoop_ode: MAX_ITER (%d) reached' ...
                 ' --> simulation may be incomplete.'], cfg.MAX_ITER);
    end

    % Close the final interval, running to the last integrated sample.
    if ~isempty(t_all)
        mode_intervals(end+1) = struct( ...
            't_start', interval_start, 't_end', t_all(end), 'mode', interval_mode);
    end

    sol.t              = t_all;
    sol.X              = X_all;
    sol.te             = te_all;
    sol.Xe             = Xe_all;
    sol.mode_intervals = mode_intervals;
    sol.ode_segments   = ode_segments;

end
