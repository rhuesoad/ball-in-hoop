function run = build_run(t_all, X_all, tau_log, u_log, scn, params, cfg)
%BUILD_RUN  Assemble a simulation result into a standard struct.
%
% The output obeys the run struct schema documented in
% results/RUN_SCHEMA.md (the single source of truth for its field names
% and units) -- read that file before changing this one.
%
% Inputs:
%   t_all   : time vector                       (N x 1) [s]
%   X_all   : physical state trajectory          (N x 6)
%   tau_log : applied (saturated) motor torque   (N x 1) [N.m]
%   u_log   : controller acceleration setpoint   (N x 1) [rad/s^2],
%             ignored (may be all-zero) for an open-loop scn
%   scn     : scenario struct (docs/MODEL.md sec. 9); scn.reference, if
%             present, must already carry the .psi_ref/.psi_dot_ref
%             handles actually used for this run (build_controller.m
%             resolves scn.reference itself for TVLQR scenarios, since
%             they don't declare one directly -- see sec. 9.4)
%   params  : physical parameters used for the simulation
%   cfg     : numerical settings used for the simulation
%
% Output:
%   run : struct obeying results/RUN_SCHEMA.md

    n = length(t_all);
    has_reference = isfield(scn, 'reference') && ~isempty(scn.reference);

    if has_reference
        ref_psi_rad       = zeros(n, 1);
        ref_psi_dot_rad_s = zeros(n, 1);
        for i = 1:n
            ref_psi_rad(i)       = scn.reference.psi_ref(t_all(i));
            ref_psi_dot_rad_s(i) = scn.reference.psi_dot_ref(t_all(i));
        end
        error_rad = ref_psi_rad - X_all(:, 5);
    else
        % Open-loop scn: there is no reference to track against.
        ref_psi_rad       = nan(n, 1);
        ref_psi_dot_rad_s = nan(n, 1);
        error_rad         = nan(n, 1);
    end

    % --- Assemble (results/RUN_SCHEMA.md) ---
    run = struct( ...
        'name',               scn.name, ...
        'source',             'simulation', ...   % docs/MODEL.md sec. 10.1; import_experiment.m sets 'experiment'
        't_s',                t_all, ...
        'x',                  X_all, ...
        'tau_Nm',             tau_log, ...
        'u_rad_s2',           u_log, ...
        'ref_psi_rad',        ref_psi_rad, ...
        'ref_psi_dot_rad_s',  ref_psi_dot_rad_s, ...
        'error_rad',          error_rad, ...
        'scenario',           scn, ...
        'params',             params, ...
        'cfg',                cfg, ...
        'timestamp',          datetime('now'));
end
