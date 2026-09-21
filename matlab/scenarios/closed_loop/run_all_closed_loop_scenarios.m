function summary = run_all_closed_loop_scenarios()
%RUN_ALL_CLOSED_LOOP_SCENARIOS  Runs every closed-loop scenario under this
%                                folder, T1 through T4.
%
%   run_all_closed_loop_scenarios()
%
%   Discovery is recursive: the scenarios sit in per-task subfolders
%   (T1/ stabilization, T2/ tracking, T3/ looping, T4/ flying ball),
%   including every attempt kept for each task. Each one's figure(s) go
%   to results/figures/<scn.id>.pdf.
%
%   Runnable standalone from any directory: it puts the project on the
%   MATLAB path itself, from its own location.
%
%   Note that the T3 and T4 scenarios plan a trajectory by direct
%   collocation, so a cold cache/trajectories/ makes this run long.
%
%   Inputs: none
%
%   Outputs
%   -------
%   summary : struct with .n_total, .n_ok and .failures (throws if any
%             scenario failed, after running all the others)

    here         = fileparts(mfilename('fullpath'));
    project_root = fileparts(fileparts(here));
    addpath(genpath(project_root));

    summary = run_scenario_set({here}, 'closed-loop');
end
