function opts = planner_opts_from_controller(tc)
%PLANNER_OPTS_FROM_CONTROLLER  A scenario's controller struct reduced to the
%                              fields plan_trajectory_casadi.m reads.
%
%   opts = planner_opts_from_controller(tc)
%
%   Inputs
%   ------
%   tc : a scenario's scn.controller struct
%
%   Outputs
%   -------
%   opts : struct carrying only the planner fields tc actually sets. Fields
%          tc leaves unset stay ABSENT, so the planner applies its own
%          documented defaults rather than this function inventing them.
%
%   WHY THIS EXISTS. build_controller.m has planner_options(), which does
%   the same mapping for the closed-loop path -- but it is a local function
%   and cannot be called from a study script. The offline studies therefore
%   each need their own copy, and a copy that falls out of step fails
%   SILENTLY: a field the scenario sets and the copy omits is dropped, and
%   the study then reports on a problem the scenario did not declare. That
%   is not hypothetical. analysis/studies/T3/T3_margin_sweep.m began with a
%   hard-coded four-field struct, and adding theta_dot_max to the T3e
%   scenario changed nothing at all until the omission was found -- two
%   full sweeps produced identical numbers that were mistaken for
%   robustness. One shared copy, here, so the next study inherits the fix
%   instead of repeating the bug.
%
%   NOTE: this is deliberately NOT merged with build_controller.m's
%   planner_options(). That one also writes the defaults out explicitly
%   because its output is the trajectory cache key (trajectory_cache_key.m),
%   where "absent" and "set to the default value" must hash the same. This
%   one is for studies that call the planner directly and never touch the
%   cache. Merging them would mean picking one of two incompatible
%   contracts; keeping the difference documented is the honest option, and
%   the FIELD LIST below is the part that must stay in step with it.

    planner_fields = {'K', 'Tf_min', 'Tf_max', 'u_max', ...
                      'contact_margin', 'psi_bounds', 'theta_dot_max', ...
                      'psi_dot_max', 'n_swings', 'liftoff', ...
                      'enforce_no_slip', 'constrain_midpoints', ...
                      'constrain_theta_f'};

    opts = struct();
    for i = 1:numel(planner_fields)
        fname = planner_fields{i};
        if isfield(tc, fname) && ~isempty(tc.(fname))
            opts.(fname) = tc.(fname);
        end
    end
end
