function plot_runs(filter)
%PLOT_RUNS  Compare multiple simulation runs from saved .mat files
%
% Usage:
%   plot_runs()                                 % all runs in 'runs/'
%   plot_runs(struct('scenario', 'Stab'))       % filter by scenario name
%   plot_runs(struct('ctrl_type', 'PID'))       % filter by controller
%   plot_runs(struct('files', {'a.mat','b.mat'})) % explicit list
%   plot_runs(struct('scenario','Stab','ctrl_type','PID')) % combined
%
% The 'scenario' filter is a substring match (case-insensitive).
% This makes it forgiving: 'bottom' matches 'Stabilization at bottom (ψ=0)'.
%
% Reads the run struct schema documented in results/RUN_SCHEMA.md. This
% function used to read run.metadata/.data/.metrics -- a struct shape
% build_run.m never actually produced (docs/AUDIT.md item C5); fixed
% here to read the real fields. The metrics table this file used to
% print (t_rise, overshoot, IAE, ITAE, saturation ratio) is not
% restored: nothing in the codebase has ever computed those values
% (compute_step_metrics.m exists but was never called, and IAE/ITAE/
% saturation_ratio were never implemented anywhere), so there was never
% real data behind that table. Adding that computation now would be new
% functionality, out of scope for a Phase 2 "no behaviour change"
% architecture pass -- flagged for a future item, not fabricated here.
%
%   Inputs
%   ------
%   filter : optional struct, see the usage examples above; default = struct()
%
%   Outputs: none (creates comparison figures)

    if nargin < 1
        filter = struct();
    end

    % Anchored on the project root rather than the current directory, so
    % this works whatever folder the caller happens to be sitting in.
    project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    runs_folder  = fullfile(project_root, 'results', 'runs');

    % --- Load runs according to filter ---
    runs = load_runs(runs_folder, filter);

    if isempty(runs)
        warning('No runs match the filter. Nothing to plot.');
        return;
    end

    fprintf('\n=== Comparing %d run(s) ===\n', numel(runs));
    for k = 1:numel(runs)
        fprintf('  %d: %s\n', k, runs{k}.name);
    end
    fprintf('\n');

    % --- Generate the three comparison plots ---
    plot_tracking_comparison(runs);
    plot_error_comparison(runs);
    plot_command_comparison(runs);
end


%% ====================================================================
%  RUN LOADING
%  ====================================================================

function runs = load_runs(folder, filter)
%LOAD_RUNS  Load .mat run files from folder, optionally filtered

    % --- Resolve which files to load ---
    % save_run.m files each run under results/runs/<YYYY-MM-DD>/, so both
    % branches below search the tree instead of reading one flat folder.
    if isfield(filter, 'files') && ~isempty(filter.files)
        % Explicit list of filenames, given without their day folder
        if ischar(filter.files)
            wanted = {filter.files};
        else
            wanted = filter.files;
        end
        filepaths = cell(size(wanted));
        for k = 1:numel(wanted)
            hit = dir(fullfile(folder, '**', wanted{k}));
            if isempty(hit)
                % Kept in the list so the load below reports it by name
                % rather than the run silently vanishing from the figure.
                filepaths{k} = fullfile(folder, wanted{k});
            else
                filepaths{k} = fullfile(hit(1).folder, hit(1).name);
            end
        end
    else
        % Every .mat under folder, day subfolders included
        listing = dir(fullfile(folder, '**', '*.mat'));
        filepaths = arrayfun(@(s) fullfile(s.folder, s.name), ...
                             listing, 'UniformOutput', false);
    end

    % --- Load and filter ---
    runs = {};
    for k = 1:numel(filepaths)
        try
            data = load(filepaths{k}, 'run');
            r = data.run;
        catch err
            warning('Could not load %s: %s', filepaths{k}, err.message);
            continue;
        end

        % Apply filters
        if isfield(filter, 'scenario') && ~isempty(filter.scenario)
            if isempty(strfind(lower(r.name), lower(filter.scenario))) %#ok<STREMP>
                continue;
            end
        end

        if isfield(filter, 'ctrl_type') && ~isempty(filter.ctrl_type)
            if ~strcmpi(r.scenario.ctrl_type, filter.ctrl_type)
                continue;
            end
        end

        runs{end+1} = r;  %#ok<AGROW>
    end

    % Sort by timestamp (chronological)
    if numel(runs) > 1
        timestamps = cellfun(@(r) r.timestamp, runs);
        [~, idx] = sort(timestamps);
        runs = runs(idx);
    end
end


%% ====================================================================
%  LABEL GENERATION
%  ====================================================================

function label = make_label(run, idx)
%MAKE_LABEL  Build a legend label from controller params (docs/MODEL.md sec. 9.3)
% Format: "PID #1 (Kp=2.0, Ki=0.5, Kd=0.1)"

    if ~isfield(run.scenario, 'controller') || isempty(run.scenario.controller)
        label = sprintf('open-loop #%d', idx);
        return;
    end

    cp   = run.scenario.controller;
    ctrl = cp.type;

    switch upper(ctrl)
        case 'PID'
            params_str = sprintf('Kp=%.2g, Ki=%.2g, Kd=%.2g', ...
                                 cp.Kp, cp.Ki, cp.Kd);

        case 'LQR'
            % Q is a 4x4 matrix; show diagonal + R
            q_diag = diag(cp.Q);
            params_str = sprintf('Q=diag(%.2g,%.2g,%.2g,%.2g), R=%.2g', ...
                                 q_diag(1), q_diag(2), q_diag(3), ...
                                 q_diag(4), cp.R);

        otherwise
            params_str = '';
    end

    label = sprintf('%s #%d (%s)', ctrl, idx, params_str);
end


%% ====================================================================
%  COMPARISON PLOTS
%  ====================================================================

function plot_tracking_comparison(runs)
%PLOT_TRACKING_COMPARISON  Overlay psi(t) for all runs vs reference

    figure('Name', 'Run comparison - Tracking');
    hold on; grid on;

    cmap = lines(numel(runs));

    % --- All trajectories ---
    for k = 1:numel(runs)
        r = runs{k};
        psi_deg = rad2deg(r.x(:, 5));
        plot(r.t_s, psi_deg, '-', ...
             'Color', cmap(k, :), 'LineWidth', 1.5, ...
             'DisplayName', make_label(r, k));
    end

    % --- Reference (assume same for all runs being compared) ---
    % If references differ across runs, this draws the first one's
    r1 = runs{1};
    plot(r1.t_s, rad2deg(r1.ref_psi_rad), 'k--', ...
         'LineWidth', 1.2, 'DisplayName', '\psi_{ref}');

    xlabel('Time [s]');
    ylabel('\psi [deg]');
    title('Reference tracking — all runs');
    legend('Location', 'best', 'FontSize', 8);
end


function plot_error_comparison(runs)
%PLOT_ERROR_COMPARISON  Overlay tracking error e(t) for all runs

    figure('Name', 'Run comparison - Tracking error');
    hold on; grid on;

    cmap = lines(numel(runs));

    for k = 1:numel(runs)
        r = runs{k};
        e_deg = rad2deg(r.error_rad);
        plot(r.t_s, e_deg, '-', ...
             'Color', cmap(k, :), 'LineWidth', 1.5, ...
             'DisplayName', make_label(r, k));
    end

    yline(0, 'k:', 'HandleVisibility', 'off');

    xlabel('Time [s]');
    ylabel('e = \psi_{ref} - \psi  [deg]');
    title('Tracking error — all runs');
    legend('Location', 'best', 'FontSize', 8);
end


function plot_command_comparison(runs)
%PLOT_COMMAND_COMPARISON  Overlay tau(t) with saturation limits

    figure('Name', 'Run comparison - Control effort');
    hold on; grid on;

    cmap = lines(numel(runs));

    % --- Saturation limits (assume same tau_max for all runs) ---
    tau_max = runs{1}.cfg.tau_max;

    for k = 1:numel(runs)
        r = runs{k};
        plot(r.t_s, r.tau_Nm, '-', ...
             'Color', cmap(k, :), 'LineWidth', 1.5, ...
             'DisplayName', make_label(r, k));
    end

    yline( tau_max, 'r--', 'LineWidth', 1.2, ...
           'Label', '+\tau_{max}', 'HandleVisibility', 'off');
    yline(-tau_max, 'r--', 'LineWidth', 1.2, ...
           'Label', '-\tau_{max}', 'HandleVisibility', 'off');
    yline(0, 'k:', 'HandleVisibility', 'off');

    xlabel('Time [s]');
    ylabel('\tau [N\cdotm]');
    ylim([-1.2*tau_max, 1.2*tau_max]);
    title('Control effort — all runs');
    legend('Location', 'best', 'FontSize', 8);
end
