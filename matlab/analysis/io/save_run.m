function filepath = save_run(run, runs_folder)
%SAVE_RUN  Persist a simulation run to disk for later comparison
%
% Saves the run struct (built by analysis/io/build_run.m)
% to a .mat file in runs_folder, with a timestamped, human-readable name.
%
% File naming convention, one subfolder per day of execution:
%   <YYYY-MM-DD>/dd_MM_yyyy-HH-mm-ss_<ctrl>_<scenario_short>.mat
%
% Example:
%   2026-05-04/04_05_2026-15-32-08_lqr_stabilization_at_bottom.mat
%
% The day folder makes a day's runs archivable or compressible as a unit;
% plot_runs.m walks into them rather than reading one flat folder.
%
% Inputs:
%   run         - run struct (must contain .metadata and .data)
%   runs_folder - parent folder, e.g. results/runs (created if it doesn't
%                 exist, along with the day subfolder)
%
% Outputs:
%   filepath    - full path of the saved file (for logging)

    % --- One subfolder per day of execution ---
    % results/runs/<YYYY-MM-DD>/. Keeps a day's runs together so a day can
    % be compressed or archived as a unit without touching the others.
    day = run.timestamp;
    day.Format = 'yyyy-MM-dd';
    runs_folder = fullfile(runs_folder, char(day));

    % --- Ensure the folder exists ---
    if ~exist(runs_folder, 'dir')
        mkdir(runs_folder);
    end

    % --- Build a clean filename from metadata ---
    ts = run.timestamp;
    ts.Format = 'dd_MM_yyyy-HH-mm-ss';
    timestamp = char(ts);
    if isfield(run.scenario, 'controller') && ~isempty(run.scenario.controller)
        ctrl_tag = lower(run.scenario.controller.type);
    else
        ctrl_tag = 'open_loop';
    end
    scn_tag   = sanitize_filename(run.name);
    % --- Handle name collisions (multiple runs in the same minute) ---
    base_name = sprintf('%s_%s_%s', timestamp, ctrl_tag, scn_tag);
    filename  = ensure_unique_filename(runs_folder, base_name);
    filepath  = fullfile(runs_folder, filename);
    
    

    % --- Save (the variable inside the .mat is named 'run') ---
    % ode_segments holds one ODE solution struct per mode interval, and those
    % carry function handles whose captured workspace is serialised with them:
    % keeping them turns a 1 MB run into 55 MB. Nothing reads them back since
    % compare_runs.m was retired, so they are dropped on the way to disk and
    % stay available on the in-memory struct.
    if isfield(run, 'ode_segments')
        run = rmfield(run, 'ode_segments');
    end
    save(filepath, 'run', '-v7.3');
    
    fprintf('  > Run saved: %s\n', filename);
end


%% ====================================================================
%  HELPERS
%  ====================================================================

function clean = sanitize_filename(name)
%SANITIZE_FILENAME  Convert a free-text name into a safe filename
% Removes/replaces characters that are problematic on Windows/Linux/macOS

    clean = lower(name);
    
    % Replace spaces and common separators with underscores
    clean = regexprep(clean, '[\s\-]+', '_');
    
    % Remove anything that's not alphanumeric or underscore
    clean = regexprep(clean, '[^a-z0-9_]', '');
    
    % Collapse multiple underscores
    clean = regexprep(clean, '_+', '_');
    
    % Trim leading/trailing underscores
    clean = regexprep(clean, '^_|_$', '');
    
    % Truncate if too long (Windows has historical 255-char path limits)
    if length(clean) > 60
        clean = clean(1:60);
    end
end


function filename = ensure_unique_filename(folder, base_name)
%ENSURE_UNIQUE_FILENAME  Append _2, _3, ... if a file already exists
% Prevents silent overwriting when multiple runs occur in the same minute

    filename = sprintf('%s.mat', base_name);
    counter = 2;
    
    while exist(fullfile(folder, filename), 'file')
        filename = sprintf('%s_%d.mat', base_name, counter);
        counter = counter + 1;
    end
end