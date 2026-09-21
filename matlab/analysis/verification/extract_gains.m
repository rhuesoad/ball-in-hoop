function gains = extract_gains(scenario_names, opts)
%EXTRACT_GAINS  Re-derives and tabulates the controller gains of each
%               closed-loop scenario, so a change anywhere in the model
%               can be checked against what it does to the controllers.
%
%   gains = extract_gains()
%   gains = extract_gains(scenario_names)
%   gains = extract_gains(scenario_names, opts)
%
%   Why this exists. The gains are not constants written down anywhere:
%   lqr_design.m solves a Riccati equation against linearize_system.m,
%   which reads rolling_matrices.m, which reads hoop_geometry.m and
%   ball_hoop_params.m. A change to the ball radius, the O-ring geometry
%   or ball_friction silently moves every gain in the project. This
%   script makes that movement visible on demand.
%
%   It runs each scenario to obtain its gains, because the scenario files
%   build their scn struct and hand it straight to run_scenario -- there
%   is no way to reach scn without executing them. run_scenario.m stores
%   the synthesised controller in run.ctrl_params for exactly this
%   purpose. Expect it to take as long as simulating the campaign.
%
%   Inputs
%   ------
%   scenario_names : cellstr of scenario function names (no arguments,
%                    returning a run struct). Defaults to
%                    default_scenarios() below.
%   opts           : (optional) struct, any subset of
%       .verbose     print the table                      (default true)
%       .stop_on_error  rethrow instead of recording the failure
%                                                        (default false)
%
%   Outputs
%   -------
%   gains : struct array, one entry per scenario, with
%       .name, .type            scenario name, controller type
%       .ok, .error             whether it ran, and why not
%       .mode, .psi_lin_rad     linearisation point (LQR only)
%       .K, .u_eq               constant gain and feedforward (LQR only)
%       .poles                  closed-loop eigenvalues (LQR only)
%       .K_schedule             n x 4 gain trajectory (TVLQR only)
%       .K_first, .K_last       its endpoints, for a quick eyeball
%       .extra                  raw ctrl_params for anything else (T4),
%                               whose gains are not a single matrix and
%                               must be read structurally
%
%   NOT a regression test. It reports what the gains are; it asserts
%   nothing about what they should be. Locking them to stored reference
%   values is refactorisation/compare_golden.m's job, and deliberately separate:
%   a gain that moved is not by itself a bug, it is a consequence to
%   look at.

    if nargin < 1 || isempty(scenario_names)
        scenario_names = default_scenarios();
    end
    if nargin < 2 || isempty(opts)
        opts = struct();
    end
    opts = apply_defaults(opts, struct('verbose', true, 'stop_on_error', false));

    if ischar(scenario_names) || isstring(scenario_names)
        scenario_names = cellstr(scenario_names);
    end

    gains = repmat(empty_gain(), 1, numel(scenario_names));
    for i = 1:numel(scenario_names)
        gains(i) = extract_one(scenario_names{i}, opts);
    end

    if opts.verbose
        print_table(gains);
    end
end


function names = default_scenarios()
%DEFAULT_SCENARIOS  The four closed-loop tasks, outer and inner where both exist.
%
%   T1_stabilization.m and T1_stabilization_inner_0deg.m are deliberately
%   NOT here: they start the ball at psi = 0 with reference 0, i.e. exactly
%   at the equilibrium being regulated to, so the solver takes 11 steps and
%   nothing is exercised. The _10deg variants are the same task posed so
%   that the controller has to do something.
    names = { ...
        'T1_stabilization_10deg', ...
        'T1_stabilization_inner_10deg', ...
        'T2_tracking_5deg_1hz', ...
        'T2_tracking_inner_5deg_1hz', ...
        'T3_loop_the_loop', ...
        'T4_flying_ball'};
end


function g = extract_one(name, opts)
%EXTRACT_ONE  Runs one scenario and reads its controller back out.
    g = empty_gain();
    g.name = name;

    try
        % evalc suppresses the scenario's own console output (lqr_design
        % and the trajectory planner are both chatty); the gains are read
        % from the returned struct, not from what they printed.
        run_struct = [];
        evalc(sprintf('run_struct = %s();', name));
    catch ME
        if opts.stop_on_error
            rethrow(ME);
        end
        g.error = ME.message;
        return;
    end

    if ~isfield(run_struct, 'ctrl_params') || isempty(run_struct.ctrl_params)
        g.error = 'run has no ctrl_params (open-loop scenario, or run_scenario.m predates storing them)';
        return;
    end

    cp   = run_struct.ctrl_params;
    g.ok = true;
    if isfield(run_struct, 'scenario') && isfield(run_struct.scenario, 'controller')
        g.type = run_struct.scenario.controller.type;
    end

    if isfield(cp, 'K') && isnumeric(cp.K) && isequal(size(cp.K), [1 4])
        g.K           = cp.K;
        g.u_eq        = get_field(cp, 'u_eq', NaN);
        g.mode        = get_field(cp, 'mode', '');
        g.psi_lin_rad = get_field(cp, 'psi_lin', NaN);
        g.poles       = get_field(cp, 'eig_cl', []);
    elseif isfield(cp, 'ref_traj') && isstruct(cp.ref_traj) && isfield(cp.ref_traj, 'K')
        K = cp.ref_traj.K;
        % Accept either n x 4 or 4 x n; report row-per-timestep.
        if size(K, 2) ~= 4 && size(K, 1) == 4
            K = K.';
        end
        g.K_schedule = K;
        g.K_first    = K(1, :);
        g.K_last     = K(end, :);
    else
        % T4 and anything else composite: no single gain matrix exists.
        % Hand back the struct rather than inventing a summary of it.
        g.extra = cp;
    end
end


function print_table(gains)
%PRINT_TABLE  Human-readable report, one block per scenario.
    fprintf('\n');
    fprintf('=== Controller gains, re-derived from the current model ===\n');
    fprintf('params/geometry in force: ball_hoop_params.m + hoop_geometry.m\n\n');

    for i = 1:numel(gains)
        g = gains(i);
        if ~g.ok
            fprintf('%-30s [FAILED] %s\n\n', g.name, g.error);
            continue;
        end

        fprintf('%-30s (%s)\n', g.name, g.type);

        if ~isempty(g.K)
            fprintf('  mode      : %s   psi_lin = %.4f rad (%.1f deg)\n', ...
                g.mode, g.psi_lin_rad, rad2deg(g.psi_lin_rad));
            fprintf('  K         : [%+11.4f %+11.4f %+11.4f %+11.4f]\n', g.K);
            fprintf('              (theta, theta_dot, psi, psi_dot) -> u [rad/s^2]\n');
            fprintf('  u_eq      : %+.4f rad/s^2\n', g.u_eq);
            if ~isempty(g.poles)
                fprintf('  poles     : %s\n', format_poles(g.poles));
            end

        elseif ~isempty(g.K_schedule)
            K = g.K_schedule;
            fprintf('  TIME-VARYING gain: %d timesteps x 4 states.\n', size(K, 1));
            fprintf('  K(t_0)    : [%+11.4f %+11.4f %+11.4f %+11.4f]\n', g.K_first);
            fprintf('  K(t_end)  : [%+11.4f %+11.4f %+11.4f %+11.4f]\n', g.K_last);
            fprintf('  range     : min [%+10.3f %+10.3f %+10.3f %+10.3f]\n', min(K, [], 1));
            fprintf('              max [%+10.3f %+10.3f %+10.3f %+10.3f]\n', max(K, [], 1));
            fprintf('  A single number cannot represent this controller -- porting it\n');
            fprintf('  means shipping the whole schedule plus its time base.\n');

        else
            fprintf('  Composite controller, no single gain matrix.\n');
            fprintf('  ctrl_params fields: %s\n', strjoin(fieldnames(g.extra)', ', '));
            fprintf('  Read it structurally before porting -- see the controller source.\n');
        end
        fprintf('\n');
    end
end


function s = format_poles(p)
%FORMAT_POLES  Compact complex-aware pole listing.
    parts = cell(1, numel(p));
    for i = 1:numel(p)
        if abs(imag(p(i))) < 1e-9
            parts{i} = sprintf('%.3f', real(p(i)));
        else
            parts{i} = sprintf('%.3f%+.3fi', real(p(i)), imag(p(i)));
        end
    end
    s = strjoin(parts, ', ');
end


function g = empty_gain()
%EMPTY_GAIN  Prototype entry, so the struct array has a fixed shape.
    g = struct( ...
        'name',        '', ...
        'type',        '', ...
        'ok',          false, ...
        'error',       '', ...
        'mode',        '', ...
        'psi_lin_rad', NaN, ...
        'K',           [], ...
        'u_eq',        NaN, ...
        'poles',       [], ...
        'K_schedule',  [], ...
        'K_first',     [], ...
        'K_last',      [], ...
        'extra',       struct());
end


function v = get_field(s, name, default_value)
%GET_FIELD  s.(name) if present and non-empty, else default_value.
    if isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default_value;
    end
end


function opts = apply_defaults(opts, defaults)
%APPLY_DEFAULTS  Fills unset fields, rejecting unknown ones loudly so a
%                misspelled option is not silently ignored.
    known = fieldnames(defaults);
    given = fieldnames(opts);
    unknown = setdiff(given, known);
    if ~isempty(unknown)
        error('extract_gains:unknown_option', ...
            'Unknown option(s): %s. Known: %s.', ...
            strjoin(unknown', ', '), strjoin(known', ', '));
    end
    for i = 1:numel(known)
        if ~isfield(opts, known{i})
            opts.(known{i}) = defaults.(known{i});
        end
    end
end
