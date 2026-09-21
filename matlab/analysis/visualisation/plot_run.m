function plot_run(run, opts)
%PLOT_RUN  Single entry point for visualizing a closed-loop run.
%
%   plot_run(run)
%   plot_run(run, opts)
%
%   Dispatches to the individual plot_* helpers in
%   analysis/visualisation/, each of which keeps doing one thing
%   (Phase 2 architecture target). Reads the run struct schema
%   documented in results/RUN_SCHEMA.md.
%
%   Inputs
%   ------
%   run  : struct built by build_run.m (results/RUN_SCHEMA.md)
%   opts : (optional) struct controlling which panels to show
%            .params              : physical parameters
%                                    (default: run.params)
%            .show_summary        : 4-panel tracking/error/torque/
%                                    phase-portrait overview
%                                    (default: true; this is what
%                                    plot_closed_loop_results.m shows)
%            .show_states         : full 6-state time series, one
%                                    subplot per state (default: false)
%            .show_phase_portrait : standalone colour-graded phase
%                                    portrait, in addition to the small
%                                    one already in the summary panel
%                                    (default: false)
%
%   Panels intentionally not wired in here, because the run struct
%   doesn't currently carry the data they need and inventing that data
%   is out of scope for this refactor (see results/RUN_SCHEMA.md and
%   docs/AUDIT.md item C5's resolution):
%     - plot_tracking.m / plot_error.m / plot_command.m need a
%       performance-metrics struct (rise time, overshoot, IAE, ITAE,
%       saturation ratio) that nothing in the codebase computes yet.
%     - plot_normal_force.m needs event times/states (te/Xe), which
%       ball_hoop_ode.m's sol struct has but build_run.m doesn't
%       currently pass through into the run struct.

    if nargin < 2, opts = struct(); end
    if ~isfield(opts, 'params'),              opts.params = run.params;         end
    if ~isfield(opts, 'show_summary'),        opts.show_summary = true;         end
    if ~isfield(opts, 'show_states'),         opts.show_states = false;         end
    if ~isfield(opts, 'show_phase_portrait'), opts.show_phase_portrait = false; end

    if opts.show_summary
        plot_closed_loop_results(run, opts.params);
    end

    if opts.show_states
        plot_states(run.t_s, run.x(:,1), run.x(:,2), run.x(:,3), ...
                    run.x(:,4), run.x(:,5), run.x(:,6));
    end

    if opts.show_phase_portrait
        psi_target = run.ref_psi_rad(end);
        if isfield(run.scenario, 'controller') && ~isempty(run.scenario.controller)
            ctrl_type = run.scenario.controller.type;
        else
            ctrl_type = 'open-loop';
        end
        plot_phase_portrait(run.x, psi_target, ctrl_type);
    end
end
