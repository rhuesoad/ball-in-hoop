function out_path = export_t4_for_bench(scn_fun, out_path)
    if nargin < 1 || isempty(scn_fun), scn_fun = @T4_flying_ball; end
    project_root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(project_root));

    run = scn_fun();
    params = run.params;
    cfg    = run.cfg;
    cp     = run.ctrl_params;

    rt     = cp.ref_traj;
    t      = rt.t_traj(:);
    x_ref  = rt.x_traj;
    u_ff   = rt.u_traj(:);
    K_traj = rt.K_traj;
    Tf     = rt.Tf;

    K_catch = cp.lqr.K(:).';
    if numel(K_catch) ~= 4
        error('export_t4_for_bench:gain', 'Phase-3 gain must have 4 entries.');
    end

    % --- MATLAB -> bench sign convention ---
    x_ref(:,3)  = -x_ref(:,3);
    x_ref(:,4)  = -x_ref(:,4);
    K_traj(:,3) = -K_traj(:,3);
    K_traj(:,4) = -K_traj(:,4);
    K_catch(3)  = -K_catch(3);
    K_catch(4)  = -K_catch(4);

    [~, R_out] = hoop_geometry('rolling_out',        params);
    [~, R_in]  = hoop_geometry('rolling_in_inside',  params);

    theta_park = pi - params.hole_center_angle;

    meta = struct( ...
        'scenario',          run.scenario.id, ...
        'Tf',                Tf, ...
        'max_abs_u',         max(abs(u_ff)), ...
        'max_abs_theta_dot', max(abs(x_ref(:,2))), ...
        'max_abs_psi_dot',   max(abs(x_ref(:,4))), ...
        'psi_release_deg',   rad2deg(x_ref(end,3)), ...
        'theta_dot_release', x_ref(end,2), ...
        'tau_max',           cfg.tau_max, ...
        'total_inertia',     params.total_inertia, ...
        'ball_radius',       params.ball_radius, ...
        'hole_width',        params.hole_angular_width, ...
        'hole_center',       params.hole_center_angle, ...
        'exported',          datestr(now, 'yyyy-mm-dd HH:MM:SS'));

    if nargin < 2 || isempty(out_path)
        out_path = fullfile(project_root, 'export', ...
                            sprintf('%s_plan.mat', run.scenario.id));
    end
    if ~exist(fileparts(out_path), 'dir'), mkdir(fileparts(out_path)); end

    save(out_path, 't', 'x_ref', 'u_ff', 'K_traj', 'Tf', ...
                   'K_catch', 'R_out', 'R_in', 'theta_park', 'meta', '-v7');

    fprintf('\nExported %s\n', out_path);
    fprintf('  %d knots, Tf = %.4f s\n', numel(t), Tf);
    fprintf('  release: psi = %.2f deg, theta_dot = %.2f rad/s\n', ...
            meta.psi_release_deg, meta.theta_dot_release);
    fprintf('  max |u_ff| = %.2f rad/s^2, max |theta_dot| = %.2f rad/s\n', ...
            meta.max_abs_u, meta.max_abs_theta_dot);
    fprintf('  R_out = %.2f mm, R_in = %.2f mm  (phase test radii)\n', ...
            1e3*R_out, 1e3*R_in);
    fprintf('  K_catch = [%.3f %.3f %.3f %.3f]\n', K_catch);
    fprintf('  theta_park = %.4f rad\n', theta_park);
end
