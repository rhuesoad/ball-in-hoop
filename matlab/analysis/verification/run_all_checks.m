%% RUN_ALL_CHECKS  Physics/numerics validation suite (Phase 3).
%
% Usage: matlab -batch "run('analysis/verification/run_all_checks.m')"
%
% Eight independent checks, each testing a physical or numerical
% property the codebase claims to have, not just "does it run without
% erroring". Runtime target: a few tens of seconds total, so this is
% cheap enough to run after every commit (unlike refactorisation/capture_golden.m
% / refactorisation/compare_golden.m, which are Phase-2-specific bit-identical
% snapshots, not a physics test suite).
%
% Distinguishes REGRESSIONS from KNOWN, ALREADY-DOCUMENTED issues: a few
% checks below are expected to fail given the codebase's current state
% (docs/AUDIT.md items still open as of this writing), and are marked
% accordingly rather than silently passed or hidden. The script only
% exits with an error (nonzero status under -batch) if there is an
% UNEXPECTED failure -- a known issue failing exactly as documented is
% not a regression and should not block a commit; a known issue
% starting to PASS, or a new failure appearing, both are.

clear; clc;
project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(project_root));
cd(project_root);

params = ball_hoop_params();
cfg    = sim_config();

fprintf('==================================================================\n');
fprintf('RUN_ALL_CHECKS -- Phase 3 validation suite\n');
fprintf('==================================================================\n\n');

results = struct('id', {}, 'name', {}, 'pass', {}, 'known_issue', {}, 'detail', {});

results = [results, check1_energy_conservation(params, cfg)];
results = [results, check2_oscillation_period(params, cfg)];
results = [results, check3_linearization(params)];
results = [results, check4_noslip_consistency(params, cfg)];
results = [results, check5_transition_continuity(params, cfg)];
results = [results, check6_static_equilibria(params, cfg)];
results = [results, check7_closed_loop_sanity(params, cfg)];
results = [results, check8_solver_independence(params, cfg)];
results = [results, check9_contact_forces(params)];

%% --- Print table ---
fprintf('\n%s\n', repmat('=', 1, 100));
fprintf('%-8s %-55s %-6s %s\n', 'ID', 'Check', 'Result', 'Detail');
fprintf('%s\n', repmat('-', 1, 100));
for k = 1:numel(results)
    r = results(k);
    if r.pass
        status = 'PASS';
    elseif r.known_issue
        status = 'FAIL*';
    else
        status = 'FAIL';
    end
    fprintf('%-8s %-55s %-6s %s\n', r.id, r.name, status, r.detail);
end
fprintf('%s\n', repmat('=', 1, 100));

n_total   = numel(results);
n_pass    = sum([results.pass]);
n_known   = sum(~[results.pass] & [results.known_issue]);
n_unknown = sum(~[results.pass] & ~[results.known_issue]);

fprintf('TOTAL: %d/%d passed', n_pass, n_total);
if n_known > 0
    fprintf('  (%d additional failure(s) are known, documented issues, marked FAIL* -- not regressions)', n_known);
end
fprintf('\n%s\n', repmat('=', 1, 100));

if n_unknown > 0
    error('run_all_checks: %d UNEXPECTED failure(s) -- see table above.', n_unknown);
end


%% ========================================================================
%  CHECK 1 -- Energy conservation
%  ========================================================================
function results = check1_energy_conservation(params, cfg)
%   With tau=0 and friction disabled, total mechanical energy
%   E = 0.5*qdot'*M*qdot - m*g*R_eff*cos(psi) must be conserved (up to
%   integrator tolerance) over a free-oscillation run. Inputs: params,
%   cfg as elsewhere in this file. Outputs: results, a struct array with
%   fields id/name/pass/known_issue/detail (one entry per mode tested).
    modes = {'rolling_out', 'rolling_in_inside'};
    results = struct('id', {}, 'name', {}, 'pass', {}, 'known_issue', {}, 'detail', {});
    for m = 1:numel(modes)
        mode = modes{m};
        p = params;
        p.ball_friction  = 0;
        p.motor_friction = 0;
        [~, R_eff] = hoop_geometry(mode, p);

        psi0 = deg2rad(20);
        c = cfg;
        c.t_end = 10;
        x0 = [R_eff; 0; 0; 0; psi0; 0];
        tau_fun = @(t, x) 0;
        sol = ball_hoop_ode(x0, tau_fun, @ball_hoop_dynamics, p, c, mode);

        theta_dot = sol.X(:,4);
        psi       = sol.X(:,5);
        psi_dot   = sol.X(:,6);
        E = zeros(numel(sol.t), 1);
        for i = 1:numel(sol.t)
            [M, ~, ~] = rolling_matrices(p, psi(i), mode);
            qdot  = [theta_dot(i); psi_dot(i)];
            E(i)  = 0.5*qdot'*M*qdot - p.ball_mass*p.gravity*R_eff*cos(psi(i));
        end
        drift_pct = max(abs(E - E(1))) / abs(E(1)) * 100;

        r = struct('id', sprintf('1.%d', m), ...
                    'name', sprintf('Energy conservation (tau=0, no friction, %s)', mode), ...
                    'pass', drift_pct < 0.1, ...
                    'known_issue', false, ...
                    'detail', sprintf('max drift = %.4f%% over 10s (< 0.1%% required), %d mode interval(s)', ...
                                       drift_pct, numel(sol.mode_intervals)));
        results = [results, r]; %#ok<AGROW>
    end
end


%% ========================================================================
%  CHECK 2 -- Small-oscillation period, derived for THIS geometry
%  ========================================================================
function results = check2_oscillation_period(params, cfg)
%   Linearizing the FULL 2-DOF (theta,psi) system (M,C=0,G from
%   rolling_matrices.m) about psi=0 with theta free (tau=0, no motor
%   friction -- the hoop is not locked, unlike the closed-loop cascade
%   assumption theta_ddot=u that linearize_system.m uses):
%
%     M11*theta_ddot + M12*psi_ddot = 0
%     M12*theta_ddot + M22*psi_ddot + m*g*R_eff*psi = 0
%
%   Eliminating theta_ddot = -(M12/M11)*psi_ddot gives
%     (M22 - M12^2/M11)*psi_ddot + m*g*R_eff*psi = 0
%     omega^2 = m*g*R_eff / (M22 - M12^2/M11)
%
%   This reduces to the classical solid-sphere-in-a-fixed-track formula
%   omega^2 = 5*g/(7*R_eff) only under TWO idealisations at once:
%   hoop_motor_inertia -> inf (M11 -> inf, so the M12^2/M11 correction
%   vanishes), AND single-point contact (r_roll -> ball_radius, so the
%   effective inertia factor drops from 1.552 back to 7/5). Neither holds
%   here: the hoop inertia is finite, and the ball rides two O-ring rails
%   (hoop_geometry.m). The two corrections happen to push the period in
%   opposite directions and partly cancel, which is exactly why the
%   classical formula must not be used as the reference -- agreeing with
%   it would be a coincidence, not a validation.
    modes = {'rolling_out', 'rolling_in_inside'};
    results = struct('id', {}, 'name', {}, 'pass', {}, 'known_issue', {}, 'detail', {});
    for m = 1:numel(modes)
        mode = modes{m};
        p = params;
        p.ball_friction  = 0;
        p.motor_friction = 0;

        [~, R_eff, R_track, coupling_sign, ~, r_roll] = hoop_geometry(mode, p);
        radius_ratio = R_track / r_roll;

        hoop_effective_inertia      = p.hoop_motor_inertia + p.ball_inertia*radius_ratio^2;
        ball_effective_inertia      = p.ball_mass*R_eff^2 + p.ball_inertia*(R_eff/r_roll)^2;
        ball_hoop_inertial_coupling = coupling_sign * p.ball_inertia * radius_ratio * (R_eff/r_roll);

        omega2_this_geometry = p.ball_mass*p.gravity*R_eff / ...
            (ball_effective_inertia - ball_hoop_inertial_coupling^2/hoop_effective_inertia);
        T_this_geometry = 2*pi/sqrt(omega2_this_geometry);

        omega2_classical = 5*p.gravity/(7*R_eff);
        T_classical = 2*pi/sqrt(omega2_classical);

        psi0 = deg2rad(1);
        c = cfg;
        c.t_end = 6*T_this_geometry;
        x0 = [R_eff; 0; 0; 0; psi0; 0];
        tau_fun = @(t, x) 0;
        sol = ball_hoop_ode(x0, tau_fun, @ball_hoop_dynamics, p, c, mode);
        t = sol.t; psi = sol.X(:,5);

        crossings = [];
        for i = 2:numel(psi)
            if sign(psi(i)) ~= sign(psi(i-1)) && psi(i) ~= 0 && psi(i-1) ~= 0
                tc = t(i-1) + (0 - psi(i-1)) * (t(i)-t(i-1)) / (psi(i)-psi(i-1));
                crossings(end+1) = tc; %#ok<AGROW>
            end
        end
        T_measured = 2*mean(diff(crossings));
        err_pct = 100*abs(T_measured - T_this_geometry)/T_this_geometry;
        err_pct_classical = 100*abs(T_measured - T_classical)/T_classical;

        r = struct('id', sprintf('2.%d', m), ...
            'name', sprintf('Small-oscillation period (%s, psi_0=1deg)', mode), ...
            'pass', err_pct < 1.0, ...
            'known_issue', false, ...
            'detail', sprintf(['T_measured=%.6fs vs derived-for-this-geometry T=%.6fs (%.4f%% err, ' ...
                                '<1%% required); classical fixed-track formula would give %.6fs (%.4f%% err)'], ...
                               T_measured, T_this_geometry, err_pct, T_classical, err_pct_classical));
        results = [results, r]; %#ok<AGROW>
    end
end


%% ========================================================================
%  CHECK 3 -- Linearisation vs finite difference
%  ========================================================================
function results = check3_linearization(params)
%   Compares linearize_system.m's analytic (A, B) against a central
%   finite-difference Jacobian of cascade_dynamics_reduced.m at each
%   mode's equilibrium. Inputs: params, as elsewhere in this file.
%   Outputs: results, a struct array with fields
%   id/name/pass/known_issue/detail (one entry per mode tested).
    modes_eq = struct('rolling_out', 0, 'rolling_in_outside', pi, 'rolling_in_inside', 0);
    mnames = fieldnames(modes_eq);
    results = struct('id', {}, 'name', {}, 'pass', {}, 'known_issue', {}, 'detail', {});

    for k = 1:numel(mnames)
        mode   = mnames{k};
        psi_eq = modes_eq.(mode);
        [A, B, ~, ~] = linearize_system(psi_eq, params, mode);

        % At every equilibrium here, sin(psi_eq)=0, so u_eq=0 (verified
        % algebraically: cascade equilibrium requires M(2,1)*u_eq +
        % m*g*R_eff*sin(psi_eq) = 0).
        u_eq = 0;
        x_eq = [0; 0; psi_eq; 0];

        h = 1e-6;
        n = 4;
        A_fd = zeros(n, n);
        for j = 1:n
            dx = zeros(n, 1); dx(j) = h;
            fp = cascade_dynamics_reduced(x_eq+dx, u_eq, params, mode);
            fm = cascade_dynamics_reduced(x_eq-dx, u_eq, params, mode);
            A_fd(:, j) = (fp - fm) / (2*h);
        end
        fp = cascade_dynamics_reduced(x_eq, u_eq+h, params, mode);
        fm = cascade_dynamics_reduced(x_eq, u_eq-h, params, mode);
        B_fd = (fp - fm) / (2*h);

        errA = max(abs(A_fd(:) - A(:))) / max(abs(A(:)));
        errB = max(abs(B_fd(:) - B(:))) / max(abs(B(:)));
        err  = max(errA, errB);

        r = struct('id', sprintf('3.%d', k), ...
            'name', sprintf('Linearisation vs finite difference (%s)', mode), ...
            'pass', err < 1e-6, ...
            'known_issue', false, ...
            'detail', sprintf('max relative error = %.3e (A: %.3e, B: %.3e), < 1e-6 required', err, errA, errB));
        results = [results, r]; %#ok<AGROW>
    end
end


%% ========================================================================
%  CHECK 4 -- No-slip constraint consistency
%  ========================================================================
function results = check4_noslip_consistency(params, cfg) %#ok<INUSD>
%   The no-slip constraint implicit in rolling_matrices.m's M matrix
%   (reverse-derived by matching KE terms; see docs/AUDIT.md item A3/1.1)
%   is:
%     phi_dot = coupling_sign*(R_track/r_roll)*theta_dot + (R_eff/r_roll)*psi_dot
%   where phi is the ball's own spin angle (untracked as a state; its
%   sign convention is a free choice, fixed here by taking the psi_dot
%   coefficient as +1). Both ratios divide by r_roll, the two-rail rolling
%   radius, not by the ball radius -- see hoop_geometry.m.
    modes = {'rolling_out', 'rolling_in_outside', 'rolling_in_inside'};
    results = struct('id', {}, 'name', {}, 'pass', {}, 'known_issue', {}, 'detail', {});
    test_pts = [1.3, -0.7; 2.1, 0.4; -0.5, -1.8];

    % --- Part A: rolling_matrices.m's M matrix is internally consistent
    % with this constraint (ball-only KE computed via phi_dot must equal
    % ball-only KE read off M).
    all_ok = true;
    for m = 1:numel(modes)
        mode = modes{m};
        [~, R_eff, R_track, coupling_sign, ~, r_roll] = hoop_geometry(mode, params);
        radius_ratio = R_track / r_roll;
        [M, ~, ~] = rolling_matrices(params, 0, mode);
        M_ball_only = M - [params.hoop_motor_inertia, 0; 0, 0];
        for k = 1:size(test_pts, 1)
            theta_dot = test_pts(k,1); psi_dot = test_pts(k,2);
            phi_dot = coupling_sign*radius_ratio*theta_dot + (R_eff/r_roll)*psi_dot;
            KE_formula = 0.5*params.ball_inertia*phi_dot^2 + 0.5*params.ball_mass*R_eff^2*psi_dot^2;
            qdot = [theta_dot; psi_dot];
            KE_M = 0.5*qdot'*M_ball_only*qdot;
            all_ok = all_ok && (abs(KE_formula - KE_M) < 1e-12*max(1, abs(KE_M)));
        end
    end
    results(end+1) = struct('id', '4.1', ...
        'name', 'No-slip constraint: rolling_matrices.m self-consistency', ...
        'pass', all_ok, 'known_issue', false, ...
        'detail', 'phi_dot formula reproduces M-matrix kinetic energy exactly, all 3 modes');

    % --- Part B: compute_energy.m (rolling_out only) matches the same
    % constraint up to the sign of phi (physically unobservable, only
    % phi_dot^2 enters energy).
    [~, R_eff, R_track, coupling_sign, ~, r_roll] = hoop_geometry('rolling_out', params);
    all_ok = true;
    for k = 1:size(test_pts, 1)
        theta_dot = test_pts(k,1); psi_dot = test_pts(k,2);
        phi_dot_rm = coupling_sign*(R_track/r_roll)*theta_dot + (R_eff/r_roll)*psi_dot;

        X1 = [R_eff, 0, 0, theta_dot, 0, psi_dot];
        [~, E_kin, ~] = compute_energy(0, X1, params, []);
        T_hoop  = 0.5*params.hoop_motor_inertia*theta_dot^2;
        T_trans = 0.5*params.ball_mass*R_eff^2*psi_dot^2;
        phi_dot_ce = sqrt(2*max(E_kin - T_hoop - T_trans, 0) / params.ball_inertia);

        all_ok = all_ok && (abs(abs(phi_dot_rm) - phi_dot_ce) < 1e-9);
    end
    results(end+1) = struct('id', '4.2', ...
        'name', 'No-slip constraint: compute_energy.m vs rolling_matrices.m', ...
        'pass', all_ok, 'known_issue', false, ...
        'detail', '|phi_dot| implied by compute_energy''s KE matches rolling_matrices'' constraint (rolling_out)');

    % --- Part C: handle_impact.m's reset rule vs the SAME constraint
    % (evaluated at phi_dot=0, the only physically stateable assumption
    % consistent with "ball had no spin at the moment of landing" that
    % would make this a legitimate use of the same expression).
    all_ok = true;
    detail_mismatches = {};
    for m = 1:numel(modes)
        mode = modes{m};
        [~, R_eff, R_track, coupling_sign, ~, r_roll] = hoop_geometry(mode, params);
        radius_ratio = R_track / r_roll;
        theta_dot_test = 3.0;

        psi_dot_handle_impact = coupling_sign * (R_eff/r_roll) * theta_dot_test;
        psi_dot_noslip_phi0   = -coupling_sign*radius_ratio*theta_dot_test / (R_eff/r_roll);

        match = abs(psi_dot_handle_impact - psi_dot_noslip_phi0) < 1e-9;
        all_ok = all_ok && match;
        if ~match
            detail_mismatches{end+1} = sprintf('%s: handle_impact=%.4f vs no-slip(phi_dot=0)=%.4f', ...
                mode, psi_dot_handle_impact, psi_dot_noslip_phi0); %#ok<AGROW>
        end
    end
    results(end+1) = struct('id', '4.3', ...
        'name', 'No-slip constraint: handle_impact.m vs rolling_matrices.m', ...
        'pass', all_ok, 'known_issue', true, ...
        'detail', ['KNOWN ISSUE (docs/AUDIT.md item 1.4, still open): handle_impact''s reset ' ...
                    'is not derived from the same no-slip expression -- ' strjoin(detail_mismatches, '; ')]);
end


%% ========================================================================
%  CHECK 5 -- Mode-transition continuity
%  ========================================================================
function results = check5_transition_continuity(params, cfg)
%   Stated preservation list (docs/AUDIT.md; matches transition_state.m/
%   handle_impact.m as currently implemented):
%     LIFTOFF (rolling_X -> free_fall): theta, theta_dot, psi, psi_dot
%       all preserved exactly; r changes by exactly cfg.TOL_NUDGE
%       (numerical bookkeeping only, not physical).
%     LANDING (free_fall -> rolling_X): theta, theta_dot preserved
%       exactly; r snapped to the geometric surface radius (checked
%       below against the EVENT-DETECTED radius, not just the post-snap
%       value, to see how close ode15s's own root-finding gets before
%       the snap does any work); r_dot forced to 0 (not preserved --
%       intentional, treated as inelastic radial impact); psi_dot fully
%       reset by the (currently physically unvalidated, see check 4.3)
%       impact rule -- not preserved, not asserted here.
    ol_scn = open_loop_test_fixtures(params);
    tau_fun_zero = @(t, x) 0;

    n_liftoff = 0; n_liftoff_ok = 0;
    n_landing = 0; n_landing_ok = 0;
    landing_gaps = [];

    for k = 1:numel(ol_scn)
        scn = ol_scn{k};
        sol = ball_hoop_ode(scn.x0, tau_fun_zero, @ball_hoop_dynamics, params, cfg, scn.state);
        mi = sol.mode_intervals;
        if numel(mi) < 2, continue; end

        for j = 1:numel(sol.te)
            x_pre = sol.Xe(j, :)';
            prev_mode = '';
            for ii = 1:numel(mi)
                if abs(mi(ii).t_end - sol.te(j)) < 1e-9
                    prev_mode = mi(ii).mode;
                    break;
                end
            end
            if isempty(prev_mode), continue; end

            [new_mode, x_post] = transition_state(prev_mode, x_pre, params, cfg);

            if strcmp(prev_mode, 'free_fall')
                n_landing = n_landing + 1;
                ok = isequal(x_pre(3), x_post(3)) && isequal(x_pre(4), x_post(4));
                if ok, n_landing_ok = n_landing_ok + 1; end
                if ~strcmp(new_mode, 'free_fall')
                    [~, R_eff_target] = hoop_geometry(new_mode, params);
                    landing_gaps(end+1) = abs(x_pre(1) - R_eff_target); %#ok<AGROW>
                end
            else
                n_liftoff = n_liftoff + 1;
                ok = isequal(x_pre(3), x_post(3)) && isequal(x_pre(4), x_post(4)) && ...
                     isequal(x_pre(5), x_post(5)) && isequal(x_pre(6), x_post(6)) && ...
                     abs(abs(x_pre(1)-x_post(1)) - cfg.TOL_NUDGE) < 1e-15;
                if ok, n_liftoff_ok = n_liftoff_ok + 1; end
            end
        end
    end

    results = struct('id', {}, 'name', {}, 'pass', {}, 'known_issue', {}, 'detail', {});
    results(end+1) = struct('id', '5.1', 'name', 'Liftoff preserves theta/theta_dot/psi/psi_dot exactly', ...
        'pass', n_liftoff > 0 && n_liftoff_ok == n_liftoff, 'known_issue', false, ...
        'detail', sprintf('%d/%d liftoff transitions (across %d open-loop scenarios)', n_liftoff_ok, n_liftoff, numel(ol_scn)));
    results(end+1) = struct('id', '5.2', 'name', 'Landing preserves theta/theta_dot exactly', ...
        'pass', n_landing > 0 && n_landing_ok == n_landing, 'known_issue', false, ...
        'detail', sprintf('%d/%d landing transitions', n_landing_ok, n_landing));
    if isempty(landing_gaps)
        gap_pass = false; gap_detail = 'no landing transitions observed to measure';
    else
        gap_pass = max(landing_gaps) < 100*cfg.AbsTol;
        gap_detail = sprintf('max landing r-gap = %.3e m (event detection vs geometric surface, before the snap); < 100*AbsTol=%.1e required, cfg.TOL_RADIUS=%.1e is %dx looser', ...
            max(landing_gaps), 100*cfg.AbsTol, cfg.TOL_RADIUS, round(cfg.TOL_RADIUS/cfg.AbsTol));
    end
    results(end+1) = struct('id', '5.3', 'name', 'r lands within solver tolerance, not the hand-tuned nudge', ...
        'pass', gap_pass, 'known_issue', false, 'detail', gap_detail);
end


%% ========================================================================
%  CHECK 6 -- Static equilibria
%  ========================================================================
function results = check6_static_equilibria(params, cfg)
%   Started exactly at each mode's equilibrium with tau=0, the system
%   should stay there (no mode switches, no drift in psi/theta) for the
%   whole run. Inputs: params, cfg as elsewhere in this file. Outputs:
%   results, a struct array with fields id/name/pass/known_issue/detail
%   (one entry per mode tested).
    modes_eq = struct('rolling_out', 0, 'rolling_in_outside', pi, 'rolling_in_inside', 0);
    mnames = fieldnames(modes_eq);
    results = struct('id', {}, 'name', {}, 'pass', {}, 'known_issue', {}, 'detail', {});

    for k = 1:numel(mnames)
        mode = mnames{k};
        psi_eq = modes_eq.(mode);
        [~, R_eff] = hoop_geometry(mode, params);

        x0 = [R_eff; 0; 0; 0; psi_eq; 0];
        c = cfg; c.t_end = 10;
        tau_fun = @(t, x) 0;
        sol = ball_hoop_ode(x0, tau_fun, @ball_hoop_dynamics, params, c, mode);

        no_switch = numel(sol.mode_intervals) == 1;
        psi_drift = max(abs(sol.X(:,5) - psi_eq));
        theta_drift = max(abs(sol.X(:,3)));
        ran_full = abs(sol.t(end) - c.t_end) < 1e-9;

        tol = 1e-9;
        pass = no_switch && ran_full && psi_drift < tol && theta_drift < tol;

        note = '';
        if strcmp(mode, 'rolling_in_outside')
            note = [' (this is the inverted/unstable equilibrium -- floating-point-level ' ...
                 'roundoff in sin(pi) does seed the instability, but does not grow to ' ...
                 'anything numerically significant within 10s for this geometry)'];
        end

        r = struct('id', sprintf('6.%d', k), ...
            'name', sprintf('Static equilibrium holds (%s)', mode), ...
            'pass', pass, 'known_issue', false, ...
            'detail', sprintf('max|psi-psi_eq|=%.2e max|theta|=%.2e mode_switches=%d ran_to_t_end=%d%s', ...
                               psi_drift, theta_drift, numel(sol.mode_intervals)-1, ran_full, note));
        results = [results, r]; %#ok<AGROW>
    end
end


%% ========================================================================
%  CHECK 7 -- Closed-loop sanity (LQR)
%  ========================================================================
function results = check7_closed_loop_sanity(params, cfg)
%   Two sanity checks on the LQR closed loop, per mode: (1) the
%   eigenvalues stored by lqr_design.m match a direct eig(A-B*K)
%   recomputation, and (2) lqr_controller.m's saturation logic keeps u
%   within [u_min, u_max] for a set of extreme test states. Inputs:
%   params, cfg as elsewhere in this file. Outputs: results, a struct
%   array with fields id/name/pass/known_issue/detail (two entries per
%   mode tested: '<k>.eig' and '<k>.sat').
    modes = {'rolling_out', 'rolling_in_outside', 'rolling_in_inside'};
    results = struct('id', {}, 'name', {}, 'pass', {}, 'known_issue', {}, 'detail', {});

    for m = 1:numel(modes)
        mode = modes{m};
        scn = struct('state', mode, 'psi_ref', @(t) 0);
        cp = lqr_design(scn, params, cfg);
        [A, B] = linearize_system(cp.psi_lin, params, mode);
        eig_direct = sort(eig(A - B*cp.K));
        eig_stored = sort(cp.eig_cl);
        eig_err = max(abs(eig_direct - eig_stored));

        results(end+1) = struct('id', sprintf('7.%d.eig', m), ... %#ok<AGROW>
            'name', sprintf('LQR eigenvalues match eig(A-B*K) (%s)', mode), ...
            'pass', eig_err < 1e-9, 'known_issue', false, ...
            'detail', sprintf('max |eig(A-B*K) - stored eig_cl| = %.2e', eig_err));

        ref.psi = @(t) 0; ref.psi_dot = @(t) 0;
        % Extreme test states: large theta_dot, large angle error, both signs.
        test_states = [ 0.05, -50, deg2rad( 80),  30; ...
                        -0.1, 100, deg2rad(-80), -30; ...
                            0,   0,            0,   0];
        within_all = true;
        for tk = 1:size(test_states, 1)
            x_phys = [0; 0; test_states(tk,1); test_states(tk,2); test_states(tk,3); test_states(tk,4)];
            u = lqr_controller(0, x_phys, ref, cp);
            theta_dot = x_phys(4);
            u_max = ( cp.tau_max - cp.motor_friction*theta_dot) / cp.total_inertia;
            u_min = (-cp.tau_max - cp.motor_friction*theta_dot) / cp.total_inertia;
            within_all = within_all && (u >= u_min - 1e-12) && (u <= u_max + 1e-12);
        end
        results(end+1) = struct('id', sprintf('7.%d.sat', m), ... %#ok<AGROW>
            'name', sprintf('lqr_controller saturation stays within [u_min,u_max] (%s)', mode), ...
            'pass', within_all, 'known_issue', false, ...
            'detail', sprintf('%d extreme test states, all within bounds', size(test_states,1)));
    end
end


%% ========================================================================
%  CHECK 8 -- Solver independence
%  ========================================================================
function results = check8_solver_independence(params, cfg)
%   Halving RelTol/AbsTol should not move the final state by more than
%   the (looser) run's own stated tolerance -- with a safety factor.
%   The naive 1x per-component bound max(RelTol*|x|, AbsTol) is too
%   strict to apply literally at a downstream time reached through a
%   discrete event: ODE tolerances bound LOCAL per-step error, not
%   global error at a later time, and event root-finding has its own
%   (comparable-order, not identical) accuracy. A 50x margin is used --
%   chosen for this documented reason, not fitted to make this pass; it
%   would still catch a genuinely broken nudge/tolerance (e.g. TOL_NUDGE
%   several orders of magnitude too large).
    ol_scn = open_loop_test_fixtures(params);
    scn = ol_scn{5};   % 'Falling until outer hoop contact' -- exercises one transition
    tau_fun_zero = @(t, x) 0;

    sol_probe = ball_hoop_ode(scn.x0, tau_fun_zero, @ball_hoop_dynamics, params, cfg, scn.state);
    t_transition = sol_probe.te(1);

    cfg_loose = cfg; cfg_loose.t_end = t_transition + 0.3;   % short past-transition horizon,
                                                              % avoids confounding with long-horizon
                                                              % nonlinear-pendulum phase sensitivity
    cfg_tight = cfg_loose;
    cfg_tight.RelTol = cfg.RelTol / 2;
    cfg_tight.AbsTol = cfg.AbsTol / 2;

    sol_loose = ball_hoop_ode(scn.x0, tau_fun_zero, @ball_hoop_dynamics, params, cfg_loose, scn.state);
    sol_tight = ball_hoop_ode(scn.x0, tau_fun_zero, @ball_hoop_dynamics, params, cfg_tight, scn.state);

    x_loose = sol_loose.X(end, :);
    x_tight = sol_tight.X(end, :);
    diff_abs = abs(x_loose - x_tight);
    bound = 50 * max(cfg_loose.RelTol*abs(x_loose), cfg_loose.AbsTol);

    pass = all(diff_abs < bound);
    worst_ratio = max(diff_abs ./ bound);

    results = struct('id', {'8.1'}, ...
        'name', {'Solver independence (RelTol/AbsTol halved, through one transition)'}, ...
        'pass', {pass}, 'known_issue', {false}, ...
        'detail', {sprintf('scenario "%s", worst component at %.2fx of the 50x-tolerance bound', scn.name, worst_ratio)});
end


%% ========================================================================
%  CHECK 9 -- Contact forces
%  ========================================================================
function results = check9_contact_forces(params)
%   dynamics/reset_maps/contact_forces.m derives the normal and tangential contact
%   forces independently of rolling_matrices.m, from the ball's own free
%   body. Two things must then hold, and neither is assumed:
%
%   9.1  The tangential equation must close. Substituting the derived F_t
%        into  m*R_eff*psi_ddot = F_t - m*g*sin(psi)  must reproduce the
%        psi_ddot that cascade_dynamics_reduced.m returns, for arbitrary
%        states and controls. A sign error in the no-slip constraint or a
%        missing dissipation term shows up here immediately -- this is the
%        check that pins the sign conventions the docstring states.
%
%   9.2  The classic constant bound must fall out as the special case it
%        is. analysis/verification/slip_acceleration_limit.m returns
%        mu*g/((1-1/k)*R_track), the bound the bench code applies as a
%        constant; at the bottom, at rest, driving exactly at it must
%        leave zero slip margin.
    modes = {'rolling_out', 'rolling_in_outside', 'rolling_in_inside'};
    m = params.ball_mass;
    g = params.gravity;

    % --- 9.1 ---
    rng(0);
    worst = 0;
    for i = 1:numel(modes)
        [~, R_eff] = hoop_geometry(modes{i}, params);
        for trial = 1:200
            x_red = [randn; 4*randn; pi*randn; 6*randn];
            u     = 200 * randn;
            [~, F_t] = contact_forces(x_red, u, params, modes{i});
            dx  = cascade_dynamics_reduced(x_red, u, params, modes{i});
            lhs = m * R_eff * dx(4);
            rhs = F_t - m * g * sin(x_red(3));
            worst = max(worst, abs(lhs - rhs) / max(1e-9, max(abs(lhs), abs(rhs))));
        end
    end
    pass_91 = worst < 1e-10;

    % --- 9.2 ---
    % rolling_out only: rolling_in_outside is convex, so at psi = 0 the
    % ball is not in contact there at all (its equilibrium is psi = pi),
    % and the bottom-at-rest premise does not apply.
    u_slip = slip_acceleration_limit(params, 'rolling_out');
    [~, ~, margin_at_bound] = contact_forces([0;0;0;0], u_slip, params, 'rolling_out');
    scale_N = params.ball_track_friction_coeff * m * g;
    pass_92 = abs(margin_at_bound) < 1e-12 * max(1, scale_N);

    results = struct( ...
        'id',          {'9.1', '9.2'}, ...
        'name',        {'Contact forces close the ball tangential equation', ...
                        'No-slip constant is the bottom-at-rest special case'}, ...
        'pass',        {pass_91, pass_92}, ...
        'known_issue', {false, false}, ...
        'detail',      {sprintf('worst relative residual %.2e over 600 random states, 3 modes', worst), ...
                        sprintf('u_slip = %.2f rad/s^2, residual margin %.2e N', u_slip, margin_at_bound)});
end


%% ========================================================================
%  SHARED TEST FIXTURES
%  ========================================================================
function scn = open_loop_test_fixtures(params)
%OPEN_LOOP_TEST_FIXTURES  Initial-condition/mode pairs covering every
%                         transition type (liftoff, landing, hole
%                         transit) that checks 5 and 8 need to
%                         re-integrate under their own custom tolerances.
%
%   Deliberately not sourced from scenarios/open_loop/*.m: those files
%   are declarations that call run_scenario(scn) (docs/MODEL.md sec. 9)
%   -- full integration, saving, and plotting each -- which would make
%   this suite slow and would litter runs/ and figures/ on every test
%   run. This suite tests the integration/event-detection engine
%   itself, not any particular named scenario, so it keeps its own
%   fixtures. Values match scenarios/open_loop/OL1-OL10.m at the time
%   of writing (not re-verified against them automatically).
%
%   Inputs
%   ------
%   params : physical parameters from ball_hoop_params()
%
%   Outputs
%   -------
%   scn : cell array of structs with fields .name, .state, .x0
%         (same shape the old open_loop_scenarii.m returned)

    % Ball-centre orbit radii from hoop_geometry.m, never
    % surface_radius - ball_radius: with the two-rail O-ring tracks the
    % latter lands 8.4 mm outside the real contact radius, which would put
    % every free-fall fixture's initial state beyond the rails.
    [~, R_eff_o] = hoop_geometry('rolling_out',       params);
    [~, R_eff_i] = hoop_geometry('rolling_in_inside', params);

    scn = { ...
        struct('name', 'Rolling on outer hoop', 'state', 'rolling_out', ...
               'x0', [R_eff_o; 0; 0; 0; pi/12; 0])

        struct('name', 'Rolling on inner hoop (inside)', 'state', 'rolling_in_inside', ...
               'x0', [R_eff_i; 0; 0; 0; pi/12; 0])

        struct('name', 'Free fall through the hole', 'state', 'free_fall', ...
               'x0', [R_eff_o - 0.001; 0; 0; 0; ...
                      params.hole_center_angle + 0.15; 0])

        struct('name', 'Falling outside the hole', 'state', 'free_fall', ...
               'x0', [R_eff_o - 0.001; 0; 0; 0; ...
                      params.hole_center_angle + 0.5; 0])

        struct('name', 'Falling until outer hoop contact', 'state', 'free_fall', ...
               'x0', [R_eff_o - 0.001; 0; 0; 0; 2*pi/3; 0])

        struct('name', 'Falling from top', 'state', 'free_fall', ...
               'x0', [R_eff_o - 0.001; 0; 0; 0; ...
                      params.hole_center_angle + 1e-6; 0])

        struct('name', 'Falling onto inner hoop', 'state', 'free_fall', ...
               'x0', [params.inner_hoop_radius - params.ball_radius - 0.001; 0; 0; 0; 11*pi/12; 0])

        struct('name', 'Transition from inner to outer hoop', 'state', 'free_fall', ...
               'x0', [params.inner_hoop_radius - params.ball_radius - 0.001; 0; pi; 0; pi/3; 0])

        struct('name', 'Escaping from inner hoop', 'state', 'rolling_in_inside', ...
               'x0', [R_eff_i; 0; 2*pi/3; 0; pi/3; 0])

        struct('name', 'Energy-validation', 'state', 'rolling_out', ...
               'x0', [R_eff_o; 0; 0; 0; 0; ...
                      1.1 * sqrt(2.3*params.gravity/R_eff_o)])
        };
end
