function validate_scenario(scn)
%VALIDATE_SCENARIO  Checks a scenario struct against docs/MODEL.md sec. 9.
%
%   validate_scenario(scn)
%
%   Fails loudly (errors, with a message naming the exact field) on any
%   missing required field or any field name not in the schema
%   (misspellings are exactly as dangerous as omissions -- a silently
%   ignored field reads to the author as "set" when it never took
%   effect). No physical quantity is given a silent default here: the
%   only optional fields are non-physical (scn.t_end, scn.id/.name
%   defaulting is not attempted at all -- both are required).
%
%   Inputs
%   ------
%   scn : scenario struct, see docs/MODEL.md sec. 9
%
%   Outputs: none (throws on failure)

    valid_modes = {'rolling_out', 'rolling_in_outside', 'rolling_in_inside', 'free_fall'};

    require_fields(scn, {'id', 'name', 'mode', 'x0'}, 'scn');
    allow_fields(scn, {'id', 'name', 'mode', 'x0', 't_end', ...
                        'controller', 'reference', 'expected_mode_sequence'}, 'scn');

    if ~ischar(scn.mode) || ~any(strcmp(scn.mode, valid_modes))
        error('validate_scenario:mode', ...
            'scn.mode must be one of {%s}, got ''%s''.', ...
            strjoin(valid_modes, ', '), toStr(scn.mode));
    end

    if ~isnumeric(scn.x0) || ~isequal(size(scn.x0), [6, 1])
        error('validate_scenario:x0', ...
            'scn.x0 must be a 6x1 numeric vector [r;r_dot;theta;theta_dot;psi;psi_dot], got a %s of size %s.', ...
            class(scn.x0), mat2str(size(scn.x0)));
    end

    if isfield(scn, 't_end') && (~isnumeric(scn.t_end) || ~isscalar(scn.t_end) || scn.t_end <= 0)
        error('validate_scenario:t_end', 'scn.t_end must be a positive scalar [s].');
    end

    if isfield(scn, 'controller') && ~isempty(scn.controller)
        validate_controller(scn.controller);
        has_controller = true;
    else
        has_controller = false;
    end

    % Both types plan their own trajectory and derive their own reference.
    is_self_referenced = has_controller && ...
        any(strcmp(scn.controller.type, {'TVLQR', 'T4'}));

    if has_controller && ~is_self_referenced
        if ~isfield(scn, 'reference') || isempty(scn.reference)
            error('validate_scenario:reference', ...
                'scn.reference is required for controller.type = ''%s'' (build it with a control/references/reference_*.m function).', ...
                scn.controller.type);
        end
        validate_reference(scn.reference);
    elseif isfield(scn, 'reference') && ~isempty(scn.reference)
        % TVLQR derives its own reference from the solved trajectory
        % (docs/MODEL.md sec. 9.4) -- a scenario-supplied one would be
        % silently ignored, which is exactly the kind of mistake this
        % function exists to catch.
        error('validate_scenario:reference', ...
            'scn.reference must not be set for controller.type = ''%s'': run_scenario.m derives it from the planned trajectory, so a scenario-supplied one would be silently ignored.', ...
            scn.controller.type);
    end

    if isfield(scn, 'expected_mode_sequence')
        if ~iscell(scn.expected_mode_sequence) || isempty(scn.expected_mode_sequence)
            error('validate_scenario:expected_mode_sequence', ...
                'scn.expected_mode_sequence must be a non-empty cell array of mode strings.');
        end
        for k = 1:numel(scn.expected_mode_sequence)
            m = scn.expected_mode_sequence{k};
            if ~ischar(m) || ~any(strcmp(m, valid_modes))
                error('validate_scenario:expected_mode_sequence', ...
                    'scn.expected_mode_sequence{%d} = ''%s'' is not a valid mode.', k, toStr(m));
            end
        end
    end
end


function validate_controller(ctrl)
%VALIDATE_CONTROLLER  Checks scn.controller against docs/MODEL.md sec. 9.3.
    require_fields(ctrl, {'type'}, 'scn.controller');

    switch ctrl.type
        case 'LQR'
            % An imposed gain (ctrl.K, for replaying a gain that actually
            % ran on the bench) and a synthesized one are alternatives, not
            % options: lqr_design.m never solves the Riccati equation when
            % K is set, so Q and R declared alongside it would describe a
            % synthesis that does not happen.
            allow_fields(ctrl, {'type', 'psi_lin', 'Q', 'R', 'K'}, 'scn.controller');
            if isfield(ctrl, 'K') && ~isempty(ctrl.K)
                require_fields(ctrl, {'psi_lin', 'K'}, 'scn.controller');
                if isfield(ctrl, 'Q') || isfield(ctrl, 'R')
                    error('validate_scenario:lqr_gain_and_weights', ...
                        ['scn.controller.K and scn.controller.Q/.R are mutually exclusive: ' ...
                         'with K set the gain is imposed, so the weights would be silently ignored.']);
                end
                if ~isnumeric(ctrl.K) || numel(ctrl.K) ~= 4
                    error('validate_scenario:lqr_gain', ...
                        'scn.controller.K must be a 4-element numeric gain [theta, theta_dot, psi, psi_dot], got a %s of size %s.', ...
                        class(ctrl.K), mat2str(size(ctrl.K)));
                end
            else
                require_fields(ctrl, {'psi_lin', 'Q', 'R'}, 'scn.controller');
            end

        case 'TVLQR'
            require_fields(ctrl, {'x0_red', 'xf_red', 'K', 'Tf_min', 'Tf_max', ...
                                   'u_max', 'Q', 'R', 'Qf'}, 'scn.controller');
            % contact_margin is allowed but not required: absent means 0,
            % plan_trajectory_casadi.m's own default.
            allow_fields(ctrl, {'type', 'x0_red', 'xf_red', 'K', 'Tf_min', 'Tf_max', ...
                                 'u_max', 'Q', 'R', 'Qf', 'contact_margin', ...
                                 'psi_bounds', 'constrain_theta_f', 'n_swings', 'liftoff', ...
                                 'theta_dot_max', 'psi_dot_max', 'enforce_no_slip', 'constrain_midpoints'}, 'scn.controller');

        case 'T4'
            % Phase 1 is a TVLQR plan (same fields), phase 3 a stationary
            % LQR on catch_mode (lqr_Q/lqr_R), phase 2 a hold on theta.
            require_fields(ctrl, {'x0_red', 'xf_red', 'K', 'Tf_min', 'Tf_max', ...
                                   'u_max', 'Q', 'R', 'Qf', ...
                                   'catch_mode', 'lqr_Q', 'lqr_R', ...
                                   'theta_hold', 'kp_hold', 'kd_hold', 'slew_boundary'}, 'scn.controller');
            allow_fields(ctrl, {'type', 'x0_red', 'xf_red', 'K', 'Tf_min', 'Tf_max', ...
                                 'u_max', 'Q', 'R', 'Qf', 'contact_margin', ...
                                 'psi_bounds', 'constrain_theta_f', 'n_swings', 'liftoff', ...
                                 'theta_dot_max', 'psi_dot_max', 'enforce_no_slip', 'constrain_midpoints', ...
                                 'catch_mode', 'lqr_Q', 'lqr_R', 'theta_dot_hold', 'catch_follow', ...
                                 'theta_hold', 'kp_hold', 'kd_hold', 'slew_boundary'}, 'scn.controller');

        otherwise
            error('validate_scenario:controller_type', ...
                'scn.controller.type must be one of {LQR, TVLQR, T4}, got ''%s''.', toStr(ctrl.type));
    end
end


function validate_reference(ref)
%VALIDATE_REFERENCE  Checks scn.reference against docs/MODEL.md sec. 9.4.
    require_fields(ref, {'psi_ref', 'psi_dot_ref'}, 'scn.reference');
    allow_fields(ref, {'psi_ref', 'psi_dot_ref', 'psi_ddot'}, 'scn.reference');

    if ~isa(ref.psi_ref, 'function_handle') || ~isa(ref.psi_dot_ref, 'function_handle')
        error('validate_scenario:reference_handles', ...
            'scn.reference.psi_ref and .psi_dot_ref must be function handles @(t) -> value.');
    end
    if isfield(ref, 'psi_ddot') && ~isempty(ref.psi_ddot) && ~isa(ref.psi_ddot, 'function_handle')
        error('validate_scenario:reference_handles', ...
            'scn.reference.psi_ddot, if set, must be a function handle @(t) -> value.');
    end
end


function require_fields(s, names, label)
%REQUIRE_FIELDS  Errors naming the first missing field, if any.
    for i = 1:numel(names)
        if ~isfield(s, names{i}) || isempty(s.(names{i}))
            error('validate_scenario:missing_field', '%s.%s is required and was not set.', label, names{i});
        end
    end
end


function allow_fields(s, names, label)
%ALLOW_FIELDS  Errors naming the first field not in the allowlist
%              (catches misspellings, which are otherwise silently ignored).
    present = fieldnames(s);
    for i = 1:numel(present)
        if ~any(strcmp(present{i}, names))
            error('validate_scenario:unknown_field', ...
                '%s.%s is not a recognized field (check docs/MODEL.md sec. 9 for the schema -- likely a misspelling).', ...
                label, present{i});
        end
    end
end


function s = toStr(v)
%TOSTR  Best-effort string conversion for an error message.
    if ischar(v)
        s = v;
    else
        s = class(v);
    end
end
