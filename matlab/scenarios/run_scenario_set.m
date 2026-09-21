function summary = run_scenario_set(folders, set_name)
%RUN_SCENARIO_SET  Runs every scenario found under a set of folders.
%
%   summary = run_scenario_set(folders, set_name)
%
%   The single engine behind run_all_scenarios.m,
%   run_all_open_loop_scenarios.m and run_all_closed_loop_scenarios.m:
%   those three differ only in which folders they hand over, so the
%   discovery/run/export loop lives here once instead of three times.
%
%   Discovery is recursive, because closed-loop scenarios sit one level
%   further down (closed_loop/T1/, T2/, ...) than open-loop ones. Files
%   whose name starts with 'run_all_' are skipped: the per-folder runners
%   live alongside the scenarios they run, so discovering them would make
%   a runner call itself.
%
%   Each scenario's figure(s) are exported to results/figures/<scn.id>.pdf
%   -- a stable name, so LaTeX \includegraphics never breaks when this is
%   re-run. A scenario that opens more than one figure gets
%   results/figures/<scn.id>_1.pdf, _2.pdf, ... in creation order.
%
%   Inputs
%   ------
%   folders  : cell array of absolute folder paths to search
%   set_name : label used in the progress report, e.g. 'open-loop'
%
%   Outputs
%   -------
%   summary : struct with .n_total, .n_ok and .failures (cell array of
%             scenario names that threw)
%
%   Throws if any scenario failed, after running all the others -- one
%   broken scenario must not hide the state of the rest.

    project_root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(project_root));

    figures_dir = fullfile(project_root, 'results', 'figures');
    if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end

    scenario_files = [];
    for i = 1:numel(folders)
        scenario_files = [scenario_files; dir(fullfile(folders{i}, '**', '*.m'))]; %#ok<AGROW>
    end

    if ~isempty(scenario_files)
        is_runner = startsWith({scenario_files.name}, 'run_all_');
        scenario_files = scenario_files(~is_runner);
    end

    n_ok = 0;
    failures = {};

    for i = 1:numel(scenario_files)
        fname = erase(scenario_files(i).name, '.m');
        close all;
        fprintf('\n=== %s ===\n', fname);
        try
            run = feval(fname);
            export_open_figures(figures_dir, run.scenario.id);
            n_ok = n_ok + 1;
        catch ME
            fprintf('  FAILED: %s\n', ME.message);
            failures{end+1} = fname; %#ok<AGROW>
        end
    end

    close all;
    fprintf('\n%d/%d %s scenarios regenerated their figure(s).\n', ...
        n_ok, numel(scenario_files), set_name);

    summary = struct('n_total', numel(scenario_files), 'n_ok', n_ok, ...
                     'failures', {failures});

    if ~isempty(failures)
        error('run_scenario_set:failures', 'Failed (%s): %s.', ...
            set_name, strjoin(failures, ', '));
    end
end


function export_open_figures(figures_dir, id)
%EXPORT_OPEN_FIGURES  Saves every currently-open figure to
%                     results/figures/<id>.pdf (single figure) or
%                     results/figures/<id>_k.pdf (k = 1, 2, ..., creation
%                     order, if the scenario opened more than one).
%
%   Inputs
%   ------
%   figures_dir : output directory
%   id          : scn.id, used as the base filename
%
%   Outputs: none (writes PDF files)

    figs = findobj('Type', 'figure');
    [~, order] = sort([figs.Number]);   % creation order, oldest first
    figs = figs(order);

    if numel(figs) == 1
        exportgraphics(figs(1), fullfile(figures_dir, [id '.pdf']));
    else
        for k = 1:numel(figs)
            exportgraphics(figs(k), fullfile(figures_dir, sprintf('%s_%d.pdf', id, k)));
        end
    end
end
