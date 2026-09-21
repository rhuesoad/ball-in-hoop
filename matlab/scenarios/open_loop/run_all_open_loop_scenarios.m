function summary = run_all_open_loop_scenarios()
%RUN_ALL_OPEN_LOOP_SCENARIOS  Runs every open-loop scenario in this folder.
%
%   run_all_open_loop_scenarios()
%
%   The passive-dynamics set (OL1-OL10): no controller, tau = 0, the ball
%   left to the hybrid dynamics alone. Each one's figure(s) go to
%   results/figures/<scn.id>.pdf.
%
%   Runnable standalone from any directory: it puts the project on the
%   MATLAB path itself, from its own location.
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

    summary = run_scenario_set({here}, 'open-loop');
end
