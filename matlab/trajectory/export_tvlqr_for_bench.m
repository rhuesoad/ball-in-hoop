function out_path = export_tvlqr_for_bench(scn_fun, out_path)
%EXPORT_TVLQR_FOR_BENCH  Writes a planned trajectory and its TVLQR gain
%                        schedule to a .npz-readable .mat for the bench.
%
%   out_path = export_tvlqr_for_bench()                 % T3_loop_the_loop_coulomb
%   out_path = export_tvlqr_for_bench(@T3_loop_the_loop, 'path/t3_plan.mat')
%
%   WHAT THE BENCH NEEDS, AND ONLY THAT. The Python side runs the same
%   control law as T1/T2 -- u = u_ff - K*(x - x_ref) -- so all it needs is
%   the four signals that law reads, sampled on a common time grid:
%       t       (K x 1)   [s]
%       x_ref   (K x 4)   [theta, theta_dot, psi, psi_dot]
%       u_ff    (K x 1)   [rad/s^2]
%       K_traj  (K x 4)   the TVLQR gain at each knot
%   Everything else in the run struct is simulation bookkeeping.
%
%   WRITTEN AS -V7 so scipy.io.loadmat reads it directly. No new
%   dependency on either side: MATLAB writes .mat, scipy reads .mat, and
%   the bench's own logs stay .npz as they are.
%
%   THE GRID IS THE PLANNER'S, NOT THE BENCH'S. K = 60 knots over ~3 s is
%   about 20 Hz, against the bench's 50 Hz loop. The Python side
%   interpolates (t3_run_v3.py, linear on x_ref/u_ff and on K), which is
%   the right place to do it: the planner's grid is where the collocation
%   is exact, and resampling it here would bake in an interpolation the
%   bench cannot then refine.
%
%   SIGN CONVENTION -- MATLAB TO BENCH. The two sides use opposite psi
%   conventions: the bench identifies B2 = -0.402 (t1_config_v3.py) where
%   this codebase's linearisation gives +0.402, the two models agreeing to
%   about 1 % on A21 and A22 and differing only in that sign. The mapping
%   is
%       psi_bench = -psi_MATLAB,      theta and u unchanged,
%   so psi, psi_dot and the two gains that multiply them flip sign, and
%   nothing else does. The control law u = u_ff - K*(x - x_ref) is then
%   invariant: K3*(-psi) with K3 flipped reproduces K3*psi exactly.
%
%   IT IS DONE HERE, ONCE. The inner-hoop gains in t1_config_v3.py carry
%   the same flip applied BY HAND ("inversion des signes de K3 et K4"),
%   and T3 was exported without it -- so the bench tracked a plan whose
%   psi ran the wrong way, and the ball left in the opposite direction
%   from the first samples (measured on the 10/08 campaign: u_ff =
%   +15 rad/s^2 with psi_ref rising while the measured psi fell). One
%   conversion point, inside the export, is the only arrangement in which
%   that cannot silently happen again.

    if nargin < 1 || isempty(scn_fun), scn_fun = @T3_loop_the_loop_coulomb; end
    project_root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(project_root));

    run = scn_fun();

    rt = run.ctrl_params.ref_traj;
    t      = rt.t_traj(:);
    x_ref  = rt.x_traj;
    u_ff   = rt.u_traj(:);
    K_traj = rt.K_traj;
    Tf     = rt.Tf;

    if size(x_ref, 2) ~= 4 || size(K_traj, 2) ~= 4
        error('export_tvlqr_for_bench:shape', ...
            'x_ref and K_traj must be K x 4; got %s and %s.', ...
            mat2str(size(x_ref)), mat2str(size(K_traj)));
    end

    % --- MATLAB -> bench sign convention (see header) ---
    % Applied BEFORE meta is built, so every number this function reports
    % is already in the bench's convention and its printout can be read
    % directly against what t3_run_v3.py --dry-run shows.
    x_ref(:,3)  = -x_ref(:,3);      % psi
    x_ref(:,4)  = -x_ref(:,4);      % psi_dot
    K_traj(:,3) = -K_traj(:,3);     % K3
    K_traj(:,4) = -K_traj(:,4);     % K4

    params = run.params;
    cfg    = run.cfg;

    % Bench-side limits the plan implies, exported so the Python config
    % cannot silently disagree with the trajectory it is tracking. These
    % are the numbers t3_config_v3.py checks itself against.
    meta = struct( ...
        'scenario',            run.scenario.id, ...
        'Tf',                  Tf, ...
        'max_abs_u',           max(abs(u_ff)), ...
        'max_abs_theta_dot',   max(abs(x_ref(:,2))), ...
        'max_abs_psi_dot',     max(abs(x_ref(:,4))), ...
        'psi_final',           x_ref(end,3), ...
        'tau_max',             cfg.tau_max, ...
        'total_inertia',       params.total_inertia, ...
        'motor_friction',      params.motor_friction, ...
        'exported',            datestr(now, 'yyyy-mm-dd HH:MM:SS'));

    if nargin < 2 || isempty(out_path)
        out_path = fullfile(project_root, 'export', ...
                            sprintf('%s_plan.mat', run.scenario.id));
    end
    if ~exist(fileparts(out_path), 'dir'), mkdir(fileparts(out_path)); end

    save(out_path, 't', 'x_ref', 'u_ff', 'K_traj', 'Tf', 'meta', '-v7');

    fprintf('\nExported %s\n', out_path);
    fprintf('  %d knots over Tf = %.4f s (%.1f Hz average)\n', ...
            numel(t), Tf, (numel(t)-1)/Tf);
    fprintf('  max |u_ff|      = %.2f rad/s^2\n', meta.max_abs_u);
    fprintf('  max |theta_dot| = %.2f rad/s   <-- bench THETA_DOT_MAX must exceed this\n', ...
            meta.max_abs_theta_dot);
    fprintf('  max |psi_dot|   = %.2f rad/s   <-- camera blur check\n', ...
            meta.max_abs_psi_dot);
    fprintf('  psi_final       = %.2f deg\n', rad2deg(meta.psi_final));
    fprintf('  |K| range       = [%.3f, %.3f]\n', ...
            min(abs(K_traj(:))), max(abs(K_traj(:))));
end
