function ctrl_params = lqr_design(scn, params, cfg)
%LQR_DESIGN  Offline LQR synthesis: linearizes the system around the
%            scenario equilibrium and solves the Riccati equation.
%
%   ctrl_params = lqr_design(scn, params, cfg)
%
%   Inputs
%   ------
%   scn    : scenario struct, must contain at least:
%              .state         : 'rolling_out' | 'rolling_in_outside' | 'rolling_in_inside'
%              .psi_ref       : reference angle handle  (used for linearization point)
%              .ctrl_params   : optional overrides for Q, R, psi_lin, and
%                               .K -- an imposed gain, which short-circuits
%                               the Riccati synthesis entirely (Q and R are
%                               then unused and returned empty)
%   params : physical parameters from ball_hoop_params()
%   cfg    : numerical settings (must contain cfg.tau_max)
%
%   Outputs
%   -------
%   ctrl_params : struct with fields used at run time by lqr_controller:
%                   .K              (1x4) LQR gain matrix
%                   .u_eq           feedforward acceleration [rad/s^2]
%                   .mode           the rolling mode (stored for diagnostics)
%                   .psi_lin        linearization angle [rad]
%                   .Q, .R          weight matrices (stored for diagnostics;
%                                   empty when .K was imposed)
%                   .eig_cl         closed-loop eigenvalues
%                   .total_inertia  motor+hoop inertia, for saturation      [kg.m^2]
%                   .motor_friction viscous friction coeff., for saturation [N.m.s/rad]
%                   .tau_max        torque saturation limit                 [N.m]

    %% --- Determine linearization angle ---
    % By default, linearize around the natural equilibrium of the mode.
    user = struct();
    if isfield(scn, 'ctrl_params') && ~isempty(scn.ctrl_params)
        user = scn.ctrl_params;
    end
    
    [~, ~, ~, ~, psi_eq_natural] = hoop_geometry(scn.state, params);

    if isfield(user, 'psi_lin')
        psi_lin = user.psi_lin;
    else
        psi_lin = psi_eq_natural;
    end
    
    %% --- Bryson's rule defaults ---
    % These are reasonable starting values; can be overridden via
    % scn.ctrl_params.Q and .R for fine-tuning.
    theta_max     = 10 * pi;          % hoop position: very lax (cyclic coord.)
    theta_dot_max = 10;               % hoop velocity: 10 rad/s before ball detaches
    psi_max       = 0.5;              % ball position error: ~30 deg
    psi_dot_max   = 5;                % ball velocity: ~5 rad/s
    u_max         = cfg.tau_max / params.total_inertia;
    
    Q_default = diag([1/theta_max^2, 1/theta_dot_max^2, ...
                      1/psi_max^2,   1/psi_dot_max^2]);
    R_default = 1 / u_max^2;
    
    Q = get_or_default(user, 'Q', Q_default);
    R = get_or_default(user, 'R', R_default);
    
    %% --- Linearize, then impose or synthesize the gain ---
    % A gain supplied as scn.ctrl_params.K is a gain that already ran
    % somewhere else -- typically transcribed back from the bench, to
    % replay a real trial in simulation. Synthesizing a fresh one from
    % (Q, R) would then simulate a controller that never flew. A and B are
    % computed either way: they are what turns an imposed gain into
    % closed-loop poles.
    [A, B, ~, ~] = linearize_system(psi_lin, params, scn.state);

    if isfield(user, 'K') && ~isempty(user.K)
        K      = user.K(:).';
        eig_cl = eig(A - B*K);
        % The weights above describe a synthesis that did not happen.
        % ctrl_params is provenance (run_scenario.m stores it on the run),
        % so reporting Bryson defaults next to an imposed gain would
        % misattribute where that gain came from.
        Q = [];
        R = [];
    else
        [K, P, eig_cl] = lqr(A, B, Q, R);    %#ok<ASGLU>
    end
    
    %% --- Feedforward acceleration u_eq ---
    % At static equilibrium: M(2,1)*u_eq + m*g*R_eff*sin(psi_lin) = 0
    [M_eq, ~, ~] = rolling_matrices(params, psi_lin, scn.state);
    [~, R_eff] = hoop_geometry(scn.state, params);
    
    if abs(M_eq(2,1)) > 1e-12
        u_eq = -params.ball_mass * params.gravity * R_eff * sin(psi_lin) ...
               / M_eq(2,1);
    else
        u_eq = 0;
        if abs(sin(psi_lin)) > 1e-6
            warning('lqr_design: M(2,1) is zero, cannot compute feedforward. Setting u_eq = 0.');
        end
    end
    
    %% --- Pack everything into ctrl_params ---
    ctrl_params.K              = K;
    ctrl_params.u_eq           = u_eq;
    ctrl_params.mode           = scn.state;
    ctrl_params.psi_lin        = psi_lin;
    ctrl_params.Q              = Q;
    ctrl_params.R              = R;
    ctrl_params.eig_cl         = eig_cl;
    
    % Motor parameters for saturation in the controller
    ctrl_params.total_inertia  = params.total_inertia;
    ctrl_params.motor_friction = params.motor_friction;
    ctrl_params.tau_max        = cfg.tau_max;
    
    % Stored for the dynamic feedforward (tracking): effective radius and
    % physical params needed to invert the ball equation online. The mode
    % (ctrl_params.mode, set above) is what feedforward.m passes to
    % rolling_matrices -- no separate surface field needed now that
    % rolling_matrices accepts mode directly (Phase 2 consolidation).
    ctrl_params.R_eff          = R_eff;
    ctrl_params.params         = params;
    
    %% --- Diagnostic output ---
    fprintf('  LQR designed at psi_lin = %.3f rad (%.1f deg) in mode %s\n', ...
            psi_lin, rad2deg(psi_lin), scn.state);
    fprintf('  K = [%+.4f  %+.4f  %+.4f  %+.4f]\n', K);
    fprintf('  u_eq = %+.4f rad/s^2\n', u_eq);
    fprintf('  Closed-loop eigenvalues:\n');
    disp(eig_cl);
end


%% ====================================================================
%  HELPERS
%% ====================================================================

function val = get_or_default(s, field, default_val)
%GET_OR_DEFAULT  Returns s.field if present and non-empty, else default_val.
    if isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = default_val;
    end
end