function [u_fun, ctrl_params, n_ctrl_states, ref, diagnostics] = build_controller(scn, params, cfg)
%BUILD_CONTROLLER  Builds the closed-loop control function for a scenario.
%
%   [u_fun, ctrl_params, n_ctrl_states, ref, diagnostics] = build_controller(scn, params, cfg)
%
%   Dispatches on scn.controller.type (docs/MODEL.md sec. 9.3). For
%   TVLQR, this is also where the reference trajectory is planned (with
%   caching, sec. 9.6) and where a non-converged plan is turned into a
%   hard error rather than silently handed to the TVLQR design step.
%
%   Inputs
%   ------
%   scn    : scenario struct with a non-empty scn.controller (see
%            docs/MODEL.md sec. 9.3); scn.reference too, unless
%            controller.type is 'TVLQR'
%   params : physical parameters from ball_hoop_params()
%   cfg    : numerical settings from sim_config()
%
%   Outputs
%   -------
%   u_fun         : function handle [u, dz_dt, info] = u_fun(t, x_phys, z)
%   ctrl_params   : struct passed to the underlying *_controller.m, incl.
%                   the saturation fields .tau_max [N.m], .total_inertia
%                   [kg.m^2], .motor_friction [N.m.s/rad]
%   n_ctrl_states : number of controller internal states (0 for all types)
%   ref           : struct with .psi_ref, .psi_dot_ref handles -- the
%                   scenario's own scn.reference for LQR, or the
%                   interpolated planned trajectory for TVLQR
%   diagnostics   : struct, populated with whichever of .eig_cl (LQR) or
%                   .solver_info (TVLQR) applies; empty struct otherwise

    project_root = fileparts(fileparts(mfilename('fullpath')));
    diagnostics  = struct();

    switch scn.controller.type
        case 'LQR'
            [u_fun, ctrl_params, n_ctrl_states, ref, diagnostics] = build_lqr(scn, params, cfg);

        case 'TVLQR'
            [u_fun, ctrl_params, n_ctrl_states, ref, diagnostics] = build_tvlqr(scn, params, cfg, project_root);

        case 'T4'
            [u_fun, ctrl_params, n_ctrl_states, ref, diagnostics] = build_t4(scn, params, cfg, project_root);

        otherwise
            % Unreachable if validate_scenario.m ran first, but
            % build_controller.m must not assume it did.
            error('build_controller:type', 'Unknown scn.controller.type ''%s''.', scn.controller.type);
    end
end


function ref_internal = to_controller_reference(ref)
%TO_CONTROLLER_REFERENCE  Translates docs/MODEL.md sec. 9.4's reference
%                         shape (.psi_ref/.psi_dot_ref, matched to
%                         results/RUN_SCHEMA.md's ref_psi_rad/
%                         ref_psi_dot_rad_s naming) into the .psi/.psi_dot
%                         shape lqr_controller.m actually expects. Kept as
%                         a translation here
%                         rather than renaming either side: the schema
%                         name matches what gets logged, the controller
%                         name predates it and several tests already
%                         depend on it.
    ref_internal.psi     = ref.psi_ref;
    ref_internal.psi_dot = ref.psi_dot_ref;
    if isfield(ref, 'psi_ddot') && ~isempty(ref.psi_ddot)
        ref_internal.psi_ddot = ref.psi_ddot;
    end
end


function [u_fun, ctrl_params, n_ctrl_states, ref, diagnostics] = build_lqr(scn, params, cfg)
%BUILD_LQR  LQR branch of build_controller.m.
%
%   lqr_design.m takes an older scenario shape (.state, .ctrl_params);
%   this adapts the new schema to it rather than changing lqr_design.m,
%   which several other callers still use as-is.
    ref = scn.reference;
    ref_internal = to_controller_reference(ref);

    lqr_scn.state = scn.mode;
    if isfield(scn.controller, 'K') && ~isempty(scn.controller.K)
        % Imposed gain, forwarded as-is: lqr_design.m skips the Riccati
        % solve and only linearizes, to report this gain's closed-loop
        % poles. validate_scenario.m has already ruled out Q/R being set
        % alongside it, so there is nothing else to pass.
        lqr_scn.ctrl_params = struct( ...
            'psi_lin', scn.controller.psi_lin, ...
            'K',       scn.controller.K);
    else
        lqr_scn.ctrl_params = struct( ...
            'psi_lin', scn.controller.psi_lin, ...   % always explicit -- validate_scenario.m requires it
            'Q',       scn.controller.Q, ...
            'R',       scn.controller.R);
    end

    ctrl_params = lqr_design(lqr_scn, params, cfg);
    n_ctrl_states = 0;
    u_fun = @(t, x_phys, z) lqr_controller(t, x_phys, ref_internal, ctrl_params);

    diagnostics.eig_cl = ctrl_params.eig_cl;
end


function [u_fun, ctrl_params, n_ctrl_states, ref, diagnostics] = build_t4(scn, params, cfg, project_root)
%BUILD_T4  Three-phase flying-ball branch of build_controller.m.
%
%   Reuses build_tvlqr for phase 1's plan (so the caching, convergence
%   check and gain design are the same code, not a parallel copy) and
%   lqr_design for phase 3's catch, then hands both to
%   t4_flying_ball_controller.m, which selects between them on the
%   measured state.
    tc = scn.controller;

    % --- Phase 1: reuse the TVLQR branch wholesale ---
    tvlqr_scn = scn;
    tvlqr_scn.controller = rmfield(tc, intersect(fieldnames(tc), {'lqr_Q', 'lqr_R', 'theta_hold', 'theta_dot_hold', 'catch_follow', 'kp_hold', ...
                                        'kd_hold', 'slew_boundary', 'catch_mode'}));
    tvlqr_scn.controller.type = 'TVLQR';
    [~, tv_ctrl, ~, tv_ref, diagnostics] = build_tvlqr(tvlqr_scn, params, cfg, project_root);
    ref_traj = tv_ctrl.ref_traj;

    % --- Phase 3: stationary LQR at the bottom of the catch surface ---
    lqr_scn.state = tc.catch_mode;
    lqr_scn.ctrl_params = struct('psi_lin', 0, 'Q', tc.lqr_Q, 'R', tc.lqr_R);
    lqr_ctrl = lqr_design(lqr_scn, params, cfg);

    [~, R_out]  = hoop_geometry(scn.mode,      params);
    [~, R_in_i] = hoop_geometry(tc.catch_mode, params);

    % Phase 1 hands x_phys straight to tvlqr_controller.m, so it needs the
    % SAME saturation fields build_tvlqr assembles above -- tau_max,
    % total_inertia, motor_friction, theta_dot_max and u_slip_max. This
    % branch used to build its own shorter list by hand, and every field
    % tvlqr_controller.m gained since then surfaced here as an
    % "Unrecognized field name" at the first closed-loop step, i.e. after
    % the plan had already been solved. Mirroring build_tvlqr's list is
    % what stops that recurring.
    ctrl_params = struct( ...
        'ref_traj',       ref_traj, ...
        'tau_max',        cfg.tau_max, ...
        'total_inertia',  params.total_inertia, ...
        'motor_friction', params.motor_friction, ...
        'theta_dot_max',  cfg.theta_dot_max, ...
        'u_slip_max',     cfg.u_slip_max, ...
        'ctrl_params',    struct('tau_max', cfg.tau_max, ...
                                 'total_inertia', params.total_inertia, ...
                                 'motor_friction', params.motor_friction, ...
                                 'theta_dot_max', cfg.theta_dot_max, ...
                                 'u_slip_max', cfg.u_slip_max, ...
                                 'ref_traj', ref_traj), ...
        'sat',         struct('tau_max', cfg.tau_max, ...
                              'total_inertia', params.total_inertia, ...
                              'motor_friction', params.motor_friction, ...
                              'theta_dot_max', cfg.theta_dot_max, ...
                              'u_slip_max', cfg.u_slip_max), ...
        'lqr',         lqr_ctrl, ...
        'ref_zero',    struct('psi', @(t) 0, 'psi_dot', @(t) 0), ...
        'theta_hold',  tc.theta_hold, ...
        'theta_dot_hold', optional(tc, 'theta_dot_hold', []), ...
        'catch_follow', optional(tc, 'catch_follow', []), ...
        'kp_hold',     tc.kp_hold, ...
        'kd_hold',     tc.kd_hold, ...
        'slew_boundary', tc.slew_boundary, ...
        'R_out',       R_out, ...
        'R_in_i',      R_in_i, ...
        ...% Hoop angle that puts the gap opposite the caught ball. The
        ...% gap sits at theta + hole_center_angle and phase 3 regulates
        ...% the ball to psi = 0, so the two are farthest apart at
        ...% theta = pi - hole_center_angle, modulo 2*pi. Passed in so
        ...% t4_flying_ball_controller.m does not need params.
        'theta_park_offset', pi - params.hole_center_angle, ...
        'contact_tol', cfg.TOL_RADIUS * 10);   % 1 mm: loose next to the 68 mm gap between the two radii

    n_ctrl_states = 0;   % all three phases are static state feedback
    u_fun = @(t, x_phys, z) t4_flying_ball_controller(t, x_phys, ctrl_params);

    % Logged reference (docs/MODEL.md sec. 9.4): the planned psi while the
    % ball is still on the outer hoop, then the catch target psi = 0 --
    % what the controller is actually aiming at in each phase, rather than
    % the TVLQR plan held past its end, which stops being the target the
    % moment the ball leaves the hoop.
    Tf = ref_traj.Tf;
    ref.psi_ref     = @(t) (t <= Tf) .* tv_ref.psi_ref(min(t, Tf));
    ref.psi_dot_ref = @(t) (t <= Tf) .* tv_ref.psi_dot_ref(min(t, Tf));
end


function plan_opts = planner_options(tc)
%PLANNER_OPTIONS  The scenario's controller struct reduced to exactly the
%                 fields plan_trajectory_casadi.m reads, with its defaults
%                 made explicit.
%
%   Defaults are written out here rather than left to the planner because
%   this struct is also the cache key (trajectory_cache_key.m): an option
%   left unset would otherwise hash as "absent" in one scenario and as its
%   value in another that set it explicitly to the same thing, and the two
%   would not share a cache entry despite being the same solve.
    plan_opts.K      = tc.K;
    plan_opts.Tf_min = tc.Tf_min;
    plan_opts.Tf_max = tc.Tf_max;
    plan_opts.u_max  = tc.u_max;

    plan_opts.contact_margin    = optional(tc, 'contact_margin', 0);
    plan_opts.psi_bounds        = optional(tc, 'psi_bounds', [-3.5*pi, 0.5*pi]);
    plan_opts.theta_dot_max     = optional(tc, 'theta_dot_max', 6*pi);
    plan_opts.psi_dot_max       = optional(tc, 'psi_dot_max', 6*pi);
    plan_opts.n_swings            = optional(tc, 'n_swings', 0);
    plan_opts.liftoff             = optional(tc, 'liftoff', []);
    plan_opts.enforce_no_slip     = optional(tc, 'enforce_no_slip', false);
    plan_opts.constrain_midpoints = optional(tc, 'constrain_midpoints', false);
    plan_opts.constrain_theta_f   = optional(tc, 'constrain_theta_f', false);
end


function v = optional(s, name, default_value)
%OPTIONAL  s.(name) if set and non-empty, else default_value.
    if isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default_value;
    end
end


function [u_fun, ctrl_params, n_ctrl_states, ref, diagnostics] = build_tvlqr(scn, params, cfg, project_root)
%BUILD_TVLQR  TVLQR branch of build_controller.m: plan (with caching),
%             verify convergence, design the time-varying gain.
    tc = scn.controller;

    % Every planner option is assembled here, once, and the same struct is
    % both hashed and handed to the solver -- see trajectory_cache_key.m
    % for the bug that motivated doing it in that order.
    plan_opts = planner_options(tc);

    cache_key  = trajectory_cache_key(tc.x0_red, tc.xf_red, scn.mode, params, plan_opts);
    cache_dir  = fullfile(project_root, 'cache', 'trajectories');
    cache_path = fullfile(cache_dir, [cache_key '.mat']);

    if exist(cache_path, 'file')
        cached = load(cache_path, 't_traj', 'x_traj', 'u_traj', 'Tf_opt', 'solver_info');
        t_traj = cached.t_traj; x_traj = cached.x_traj; u_traj = cached.u_traj;
        Tf_opt = cached.Tf_opt; solver_info = cached.solver_info;
    else
        [t_traj, x_traj, u_traj, Tf_opt, sol] = plan_trajectory_casadi( ...
            tc.x0_red, tc.xf_red, plan_opts, params, scn.mode);

        st = sol.stats();
        solver_info.exit_status             = st.return_status;
        solver_info.converged               = st.success;
        solver_info.iterations              = st.iter_count;
        solver_info.kkt_residual            = st.iterations.inf_du(end);
        solver_info.max_constraint_violation = st.iterations.inf_pr(end);

        if ~solver_info.converged
            error('build_controller:tvlqr_not_converged', ...
                'TVLQR trajectory planning did not converge (IPOPT status: %s). Refusing to hand an unconverged trajectory to tvlqr_design.m.', ...
                solver_info.exit_status);
        end

        if ~exist(cache_dir, 'dir'), mkdir(cache_dir); end
        save(cache_path, 't_traj', 'x_traj', 'u_traj', 'Tf_opt', 'solver_info');
    end

    K_traj = tvlqr_design(t_traj, x_traj, u_traj, tc.Q, tc.R, tc.Qf, params);

    ref_traj.t_traj = t_traj;
    ref_traj.x_traj = x_traj;
    ref_traj.u_traj = u_traj;
    ref_traj.K_traj = K_traj;
    ref_traj.Tf     = Tf_opt;

    ctrl_params.tau_max        = cfg.tau_max;
    ctrl_params.total_inertia  = params.total_inertia;
    ctrl_params.motor_friction = params.motor_friction;
    % AJOUTS POUR SLIP ET THETAMAX
    ctrl_params.theta_dot_max = cfg.theta_dot_max;
    ctrl_params.u_slip_max    = cfg.u_slip_max;
    % Exposed so a caller that composes this branch into a larger
    % controller (build_t4) can reach the planned trajectory without
    % digging it back out of the returned closure.
    ctrl_params.ref_traj       = ref_traj;

    n_ctrl_states = 0;
    u_fun = @(t, x_phys, z) tvlqr_controller(t, x_phys, ref_traj, ctrl_params);

    % Reference for logging (docs/MODEL.md sec. 9.4): the planned psi(t),
    % interpolated the same way tvlqr_controller.m holds the last knot
    % beyond Tf.
    ref.psi_ref     = @(t) interp1(t_traj, x_traj(:,3), min(max(t, t_traj(1)), Tf_opt), 'linear', 'extrap');
    ref.psi_dot_ref = @(t) interp1(t_traj, x_traj(:,4), min(max(t, t_traj(1)), Tf_opt), 'linear', 'extrap');

    diagnostics.solver_info = solver_info;
end
