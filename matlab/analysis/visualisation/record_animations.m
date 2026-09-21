function out = record_animations(varargin)
%RECORD_ANIMATIONS  Write one MP4 per saved run, into results/animations.
%
%   record_animations()
%   record_animations('runs_dir', d, 'out_dir', o, 'fps', 30, 'filter', 'T3')
%
%   Reads the runs save_run.m produced and plays each one through
%   animate_system.m with its video output on. The file is named from the
%   scenario id rather than the run's timestamp, so re-running the scenarios
%   overwrites the clip instead of piling up a second copy.
%
%   'filter' keeps only the runs whose file name contains the string, which
%   is how you re-record a single task without redoing all of them.

    project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));

    p = inputParser;
    p.addParameter('runs_dir', fullfile(project_root, 'results', 'runs'), @(s) ischar(s) || isstring(s));
    p.addParameter('out_dir',  fullfile(project_root, 'results', 'animations'), @(s) ischar(s) || isstring(s));
    p.addParameter('fps',      30, @isscalar);
    p.addParameter('filter',   '', @(s) ischar(s) || isstring(s));
    p.parse(varargin{:});
    opt = p.Results;

    files = dir(fullfile(char(opt.runs_dir), '**', '*.mat'));
    if ~isempty(char(opt.filter))
        files = files(contains({files.name}, char(opt.filter)));
    end
    if isempty(files)
        error('record_animations:no_runs', 'No runs under %s.', char(opt.runs_dir));
    end

    if ~exist(char(opt.out_dir), 'dir')
        mkdir(char(opt.out_dir));
    end

    fprintf('%d run(s) -> %s\n', numel(files), char(opt.out_dir));
    out = struct('name', {}, 'file', {}, 'frames', {}, 'seconds', {});

    for k = 1:numel(files)
        path = fullfile(files(k).folder, files(k).name);
        S = load(path, 'run');
        run = S.run;

        name = run.scenario.id;
        dest = fullfile(char(opt.out_dir), [name '.mp4']);

        fprintf('  [%2d/%2d] %-34s %5.2f s ... ', k, numel(files), name, run.t_s(end));
        animate_system(run.t_s, run.x, run.params, [], ...
                       'video', dest, 'fps', opt.fps);
        close(gcf);

        d = dir(dest);
        out(end+1) = struct('name', name, 'file', dest, ...
                            'frames', round(run.t_s(end) * opt.fps), ...
                            'seconds', run.t_s(end)); %#ok<AGROW>
        fprintf('      %.1f MB\n', d.bytes / 1e6);
    end

    fprintf('%d animation(s) written.\n', numel(out));
end
