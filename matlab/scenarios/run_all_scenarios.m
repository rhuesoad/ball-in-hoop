function summary = run_all_scenarios()
%RUN_ALL_SCENARIOS  Runs every scenario file and regenerates every
%                    thesis figure in one command.
%
%   run_all_scenarios()
%
%   Runs the open-loop set (scenarios/open_loop/, OL1-OL10) and the
%   closed-loop set (scenarios/closed_loop/, T1-T4) in that order, and
%   exports whatever figure(s) run_scenario.m produced to
%   results/figures/<scn.id>.pdf.
%
%   To run only one of the two sets, call
%   run_all_open_loop_scenarios.m or run_all_closed_loop_scenarios.m.
%
%   Inputs: none
%
%   Outputs
%   -------
%   summary : struct with .n_total, .n_ok and .failures (throws if any
%             scenario failed, after running all the others)

    scenarios_dir = fileparts(mfilename('fullpath'));

    summary = run_scenario_set( ...
        {fullfile(scenarios_dir, 'open_loop'), ...
         fullfile(scenarios_dir, 'closed_loop')}, 'all');
end
