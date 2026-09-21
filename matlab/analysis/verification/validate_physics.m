%% VALIDATE_PHYSICS  Regression checks accumulated during Phase 1 (docs/AUDIT.md fixes).
%
% Each cell below is a self-contained check for one AUDIT.md item, added in
% the commit that implements that item's fix. Run the whole script after
% any change to dynamics/ or control/ to catch regressions; each cell also
% prints the numbers quoted in that commit's message.
%
% Usage: matlab -batch "run('analysis/verification/validate_physics.m')"

clear; clc;
project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(project_root));
cd(project_root);   % save_run.m saves to the cwd-relative 'runs/'
params = ball_hoop_params();

n_pass = 0;
n_fail = 0;

fprintf('==================================================================\n');
fprintf('VALIDATE_PHYSICS -- Phase 1 regression suite\n');
fprintf('==================================================================\n\n');

%% ---------------------------------------------------------------------
%  A3/Phase2 -- hoop_geometry.m is the single source of the mode ->
%  surface/radius mapping. This check used to assert that it reproduced
%  the old inline formula R_contact -/+ ball_radius exactly, which was
%  right while that formula was the model. It no longer is: the ball
%  rides two O-ring rails, so R_eff comes from thesis eq. 78 and is
%  8.4 mm smaller on the outer track. The reference below is therefore
%  recomputed here from the O-ring dimensions, independently of
%  hoop_geometry.m -- a cross-check of the formula, not a regression
%  lock on a superseded number. The pre-two-rail values are printed
%  alongside so the size of the correction stays visible.
%% ---------------------------------------------------------------------
fprintf('--- A3/Phase2: hoop_geometry.m ---\n');

modes = {'rolling_out', 'rolling_in_outside', 'rolling_in_inside'};

% Independent restatement of eq. 78: Delta = sqrt((Rb+rc)^2 - (h/2)^2),
% R_eff = Rg -/+ Delta with - on the concave side, + on the convex one.
Delta_ref = sqrt((params.ball_radius + params.oring_cord_radius)^2 ...
                 - (params.oring_axial_spacing/2)^2);
expected_Reff = [ ...
    params.outer_hoop_oring_radius  - Delta_ref, ...
    params.inner_hoop_outer_radius  + Delta_ref, ...
    params.inner_hoop_inner_radius  - Delta_ref];
single_contact_Reff = [ ...
    params.outer_hoop_inner_radius - params.ball_radius, ...
    params.inner_hoop_outer_radius + params.ball_radius, ...
    params.inner_hoop_inner_radius - params.ball_radius];
old_surface = {'outer', 'inner_outside', 'inner_inside'};
old_psi_eq  = [0, pi, 0];

for k = 1:numel(modes)
    [surface, R_eff, ~, ~, psi_eq] = hoop_geometry(modes{k}, params);
    ok_Reff = abs(R_eff - expected_Reff(k)) < 1e-15;
    ok_surf = strcmp(surface, old_surface{k});
    ok_eq   = abs(psi_eq - old_psi_eq(k)) < 1e-15;
    ok = ok_Reff && ok_surf && ok_eq;
    fprintf('  %-20s R_eff=%.6f m (single-contact would give %.6f m)  surface=%-14s psi_eq=%.4f rad  [%s]\n', ...
        modes{k}, R_eff, single_contact_Reff(k), surface, psi_eq, pass_str(ok));
    [n_pass, n_fail] = tally(ok, n_pass, n_fail);
end

% Consumers wired to hoop_geometry still produce the same numbers as before.
[A_out, B_out] = linearize_system(0, params, 'rolling_out');
ok = abs(B_out(2) - 1.0) < 1e-12;   % cascade signature, unaffected by A3
fprintf('  linearize_system(''rolling_out'') runs, B(2,1)=%.4f  [%s]\n', B_out(2), pass_str(ok));
[n_pass, n_fail] = tally(ok, n_pass, n_fail);

scn.state = 'rolling_out';
scn.psi_ref = @(t) 0;
scn.ctrl_params = struct();
cfg = sim_config();
ctrl_params = lqr_design(scn, params, cfg);
ok = strcmp(ctrl_params.mode, 'rolling_out') && abs(ctrl_params.R_eff - expected_Reff(1)) < 1e-15;
fprintf('  lqr_design(''rolling_out'') runs, mode=%s R_eff=%.6f  [%s]\n', ...
    ctrl_params.mode, ctrl_params.R_eff, pass_str(ok));
[n_pass, n_fail] = tally(ok, n_pass, n_fail);

fprintf('\n');

%% ---------------------------------------------------------------------
%  C2/1.3 -- handle_impact.m now rotates the hole with the hoop, and
%  agrees with detect_events.m's wraparound-safe test
%% ---------------------------------------------------------------------
fprintf('--- C2/1.3: hole-boundary offset in handle_impact.m ---\n');

[~, R_eff_ii] = hoop_geometry('rolling_in_inside', params);

% (a) theta=0, psi=180deg (center of hole): unaffected by the fix --
%     both old and new formulas agree the ball is over the gap.
x_pre = [R_eff_ii; 0; 0; 0; deg2rad(180); 0];
[new_mode_a, ~] = handle_impact(x_pre, params, cfg);
ok_a = strcmp(new_mode_a, 'free_fall');
fprintf('  theta=0, psi=180deg (hole center):        mode=%-16s [%s] (unaffected by fix)\n', ...
    new_mode_a, pass_str(ok_a));
[n_pass, n_fail] = tally(ok_a, n_pass, n_fail);

% (b) theta=0, psi=0deg (opposite the hole): unaffected by the fix.
x_pre = [R_eff_ii; 0; 0; 0; deg2rad(0); 0];
[new_mode_b, ~] = handle_impact(x_pre, params, cfg);
ok_b = strcmp(new_mode_b, 'rolling_in_inside');
fprintf('  theta=0, psi=0deg (opposite hole):         mode=%-16s [%s] (unaffected by fix)\n', ...
    new_mode_b, pass_str(ok_b));
[n_pass, n_fail] = tally(ok_b, n_pass, n_fail);

% (c) theta=60deg, psi=250deg: the hole has rotated with the hoop to
%     [195,285]deg. The OLD formula ignored hoop_theta and tested against
%     the un-rotated window [135,225]deg -- psi=250deg is outside that,
%     so it would have (wrongly) classified this as rolling contact.
%     The FIXED formula correctly recognizes the ball is over the
%     (rotated) gap.
theta_c = deg2rad(60);
psi_c   = deg2rad(250);
x_pre   = [R_eff_ii; 0; theta_c; 0; psi_c; 0];
[new_mode_c, ~] = handle_impact(x_pre, params, cfg);

old_hole_start = mod(params.hole_center_angle - params.hole_angular_width/2, 2*pi);
old_hole_end   = mod(params.hole_center_angle + params.hole_angular_width/2, 2*pi);
old_in_hole    = (mod(psi_c,2*pi) >= old_hole_start - cfg.TOL_ANGLE) && ...
                 (mod(psi_c,2*pi) <= old_hole_end   + cfg.TOL_ANGLE);
if old_in_hole, old_mode_c = 'free_fall'; else, old_mode_c = 'rolling_in_inside'; end

ok_c = strcmp(new_mode_c, 'free_fall');
fprintf('  theta=60deg, psi=250deg (rotated hole):    BEFORE mode=%-16s AFTER mode=%-16s [%s]\n', ...
    old_mode_c, new_mode_c, pass_str(ok_c));
[n_pass, n_fail] = tally(ok_c, n_pass, n_fail);

fprintf('\n');

%% ---------------------------------------------------------------------
%  1.2 -- normal-force sign at the rolling_in_outside equilibrium
%% ---------------------------------------------------------------------
fprintf('--- 1.2: normal-force sign for convex contact (rolling_in_outside) ---\n');

[~, R_eff_io] = hoop_geometry('rolling_in_outside', params);
x_eq     = [R_eff_io; 0; 0; 0; pi; 0];   % at rest, at the natural equilibrium psi=pi
N_after  = detect_events(0, x_eq, params, cfg, 'rolling_in_outside');
N_before = params.ball_mass * (params.gravity * cos(pi) + R_eff_io * 0^2);  % the old (wrong) formula

ok = N_after > 0;
fprintf('  N at rest, psi=pi (equilibrium):  BEFORE=%+.4f N (spurious liftoff)  AFTER=%+.4f N  [%s]\n', ...
    N_before, N_after, pass_str(ok));
[n_pass, n_fail] = tally(ok, n_pass, n_fail);

fprintf('\n');

%% ---------------------------------------------------------------------
%  C4 -- plot_normal_force.m sign fix
%% ---------------------------------------------------------------------
fprintf('--- C4: plot_normal_force.m sign ---\n');

% At rest at the rolling_out equilibrium (psi=0), a concave-contact ball
% must show a positive (compressive) normal force. The formula itself
% (not the plotting call, which has no return value) is what changed;
% reproduce it here exactly as the fixed file now computes it.
[~, r_eq] = hoop_geometry('rolling_out', params);
N_after = params.ball_mass * (params.gravity * cos(0) + r_eq * 0^2);
ok = N_after > 0;
fprintf('  N at rest, psi=0 (both forms coincide when psi_dot=0): AFTER=%+.4f N  [%s]\n', N_after, pass_str(ok));
[n_pass, n_fail] = tally(ok, n_pass, n_fail);

% Case where the sign actually matters: fast rotation, so the centripetal
% term dominates and the (wrong) "-" formula would go negative/misleading
% while the correct "+" formula stays positive (both terms add).
psi_dot_fast = 5;   % rad/s
N_after2  = params.ball_mass * (params.gravity * cos(0) + r_eq * psi_dot_fast^2);
N_before2 = params.ball_mass * (params.gravity * cos(0) - r_eq * psi_dot_fast^2);
ok2 = N_after2 > N_before2;
fprintf('  N at psi=0, psi_dot=%.0f rad/s:  BEFORE(wrong "-")=%+.4f N  AFTER(fixed "+")=%+.4f N  [%s]\n', ...
    psi_dot_fast, N_before2, N_after2, pass_str(ok2));
[n_pass, n_fail] = tally(ok2, n_pass, n_fail);

% Smoke test: the plotting function itself still runs without error.
try
    te_all = []; Xe_all = [];
    t_smoke   = linspace(0, 1, 10).';
    r_smoke   = r_eq * ones(10, 1);
    psi_smoke = zeros(10, 1);
    psid_smoke = zeros(10, 1);
    plot_normal_force(params, t_smoke, r_smoke, psi_smoke, psid_smoke, te_all, Xe_all);
    close(gcf);
    ok3 = true;
catch ME
    ok3 = false;
    fprintf('  plot_normal_force smoke test threw: %s\n', ME.message);
end
fprintf('  plot_normal_force(...) runs without error  [%s]\n', pass_str(ok3));
[n_pass, n_fail] = tally(ok3, n_pass, n_fail);

fprintf('\n');

%% ---------------------------------------------------------------------
%  C7 -- cascade_dynamics_reduced.m accepts a mode argument
%% ---------------------------------------------------------------------
fprintf('--- C7: cascade_dynamics_reduced.m mode argument ---\n');

x_red = [0; 0.5; deg2rad(10); 0.2];
u_red = 1.0;

dx_3arg    = cascade_dynamics_reduced(x_red, u_red, params);          % old call signature
dx_default = cascade_dynamics_reduced(x_red, u_red, params, 'rolling_out');
ok = isequal(dx_3arg, dx_default);
fprintf('  3-arg call == 4-arg call with mode=''rolling_out'' (default):  [%s]\n', pass_str(ok));
[n_pass, n_fail] = tally(ok, n_pass, n_fail);

dx_inner = cascade_dynamics_reduced(x_red, u_red, params, 'rolling_in_inside');
ok2 = ~isequal(dx_3arg, dx_inner);   % different surface -> genuinely different dynamics
fprintf('  mode=''rolling_in_inside'' gives different psi_ddot than outer: dx_outer(4)=%.6f  dx_inner(4)=%.6f  [%s]\n', ...
    dx_3arg(4), dx_inner(4), pass_str(ok2));
[n_pass, n_fail] = tally(ok2, n_pass, n_fail);

fprintf('\n');

%% ---------------------------------------------------------------------
%  1.4 (partial) -- slip_sign now matches rolling_matrices' coupling_sign
%% ---------------------------------------------------------------------
fprintf('--- 1.4 (partial): slip_sign vs coupling_sign in handle_impact.m ---\n');

theta_dot_pre = 2.0;   % rad/s, arbitrary nonzero hoop rate at impact

% rolling_out: the slip_sign fix left this case alone (old hardcoded -1 and
% new coupling_sign agree), but the two-rail correction does not: the reset
% now divides by r_roll, not by the ball radius, so psi_dot after landing is
% ~18 % larger. The reference below is updated accordingly -- what is still
% being checked is that handle_impact uses hoop_geometry's rolling radius
% and not params.ball_radius.
[~, R_eff_o, ~, ~, ~, r_roll_o] = hoop_geometry('rolling_out', params);
x_pre   = [R_eff_o; 0; 0; theta_dot_pre; deg2rad(0); 0];
[~, x0] = handle_impact(x_pre, params, cfg);
psi_dot_expected     = -(R_eff_o/r_roll_o) * theta_dot_pre;
psi_dot_single_point = -(R_eff_o/params.ball_radius) * theta_dot_pre;
ok = abs(x0(6) - psi_dot_expected) < 1e-12;
fprintf('  rolling_out landing:         psi_dot AFTER=%+.4f (single-contact would give %+.4f)  [%s]\n', ...
    x0(6), psi_dot_single_point, pass_str(ok));
[n_pass, n_fail] = tally(ok, n_pass, n_fail);

% rolling_in_inside: sign was wrong (old slip_sign=+1, correct coupling_sign=-1).
[~, R_eff_ii] = hoop_geometry('rolling_in_inside', params);
x_pre    = [R_eff_ii; 0; 0; theta_dot_pre; deg2rad(0); 0];
[~, x0]  = handle_impact(x_pre, params, cfg);
psi_dot_before = (+1) * (R_eff_ii/params.ball_radius) * theta_dot_pre;   % old (wrong) slip_sign=+1
psi_dot_after  = x0(6);
ok = sign(psi_dot_after) ~= sign(psi_dot_before);
fprintf('  rolling_in_inside landing:   BEFORE psi_dot=%+.4f  AFTER psi_dot=%+.4f  [%s]\n', ...
    psi_dot_before, psi_dot_after, pass_str(ok));
[n_pass, n_fail] = tally(ok, n_pass, n_fail);

% rolling_in_outside: sign was wrong (old slip_sign=-1, correct coupling_sign=+1).
[~, R_eff_io] = hoop_geometry('rolling_in_outside', params);
x_pre    = [R_eff_io; 0; 0; theta_dot_pre; deg2rad(0); 0];
[~, x0]  = handle_impact(x_pre, params, cfg);
psi_dot_before = (-1) * (R_eff_io/params.ball_radius) * theta_dot_pre;   % old (wrong) slip_sign=-1
psi_dot_after  = x0(6);
ok = sign(psi_dot_after) ~= sign(psi_dot_before);
fprintf('  rolling_in_outside landing:  BEFORE psi_dot=%+.4f  AFTER psi_dot=%+.4f  [%s]\n', ...
    psi_dot_before, psi_dot_after, pass_str(ok));
[n_pass, n_fail] = tally(ok, n_pass, n_fail);

fprintf('\n');

%% ---------------------------------------------------------------------
%  1.5 -- total_inertia now equals hoop_motor_inertia exactly (user
%  decision: matches the "ideal inner loop ignores the ball" model)
%% ---------------------------------------------------------------------
fprintf('--- 1.5: total_inertia definition ---\n');

total_inertia_before = params.hoop_motor_inertia + params.ball_radius^2 * params.ball_mass;
total_inertia_after  = params.total_inertia;
ok = abs(total_inertia_after - params.hoop_motor_inertia) < 1e-15;
fprintf('  total_inertia:  BEFORE=%.6f kg.m^2  AFTER=%.6f kg.m^2  (hoop_motor_inertia=%.6f)  [%s]\n', ...
    total_inertia_before, total_inertia_after, params.hoop_motor_inertia, pass_str(ok));
[n_pass, n_fail] = tally(ok, n_pass, n_fail);

fprintf('\n');

%% ---------------------------------------------------------------------
%  1.6/C3 -- both u->tau sites now call motor_inner_loop.m
%% ---------------------------------------------------------------------
fprintf('--- 1.6/C3: u->tau routed through motor_inner_loop.m ---\n');

u_test         = 3.0;                 % rad/s^2
theta_dot_test = 1.5;                 % rad/s
x_phys_test    = [0; 0; 0; theta_dot_test; 0; 0];

tau_after = motor_inner_loop(u_test, x_phys_test, params, cfg);   % what both call sites now compute

% What augmented_dynamics computed before this commit (already using the
% post-1.5 total_inertia, since that landed in an earlier commit):
tau_augdyn_before = u_test * params.total_inertia + params.motor_friction * theta_dot_test;

% What call_controller_tau computed before this commit (hoop_motor_inertia,
% the 1.6 bug) -- note this already coincides with tau_after numerically
% now that 1.5 set total_inertia = hoop_motor_inertia; the remaining value
% of this commit is a single source of truth (C3), not a numeric change.
tau_logtau_before = u_test * params.hoop_motor_inertia + params.motor_friction * theta_dot_test;

ok1 = abs(tau_after - tau_augdyn_before) < 1e-15;
ok2 = abs(tau_after - tau_logtau_before) < 1e-15;
fprintf('  motor_inner_loop tau = %.6f N.m\n', tau_after);
fprintf('  matches augmented_dynamics'' pre-commit formula (total_inertia):    [%s]\n', pass_str(ok1));
fprintf('  matches call_controller_tau''s pre-commit formula (hoop_motor_inertia): [%s]\n', pass_str(ok2));
fprintf('  (both now match because 1.5 already set total_inertia = hoop_motor_inertia;\n');
fprintf('   this commit''s value is deduplication, not a further numeric change here.)\n');
[n_pass, n_fail] = tally(ok1, n_pass, n_fail);
[n_pass, n_fail] = tally(ok2, n_pass, n_fail);

fprintf('\n');

%% ---------------------------------------------------------------------
%  1.7/C6 -- plan_trajectory_casadi.m and tvlqr_design.m now match
%  rolling_matrices.m instead of an independent (buggy) re-derivation
%% ---------------------------------------------------------------------
fprintf('--- 1.7/C6: CasADi/TVLQR dynamics now match rolling_matrices.m ---\n');

% (a) tvlqr_design.m's linearize_cascade is now a thin wrapper around
% linearize_system (verified by inspection: it's a 3-line delegation,
% see tvlqr_design.m), so it is correct by construction as long as
% linearize_system itself is (already checked in the A3 section above).
% What's worth checking numerically is the size of the change from the
% old, independently hand-derived formula it replaced:
x_eq_test = [0; 0.3; deg2rad(5); 0.1];
[A_ref, B_ref]   = linearize_system(x_eq_test(3), params, 'rolling_out');

% Old (buggy) B(4) coefficient, for comparison: e/a with e=I_b*r*(r+1),
% a=m*R_eff^2+I_b*r^2, r=R_o/R_b -- the historical formula being replaced.
m_b = params.ball_mass; R_b = params.ball_radius; R_o = params.outer_hoop_inner_radius;
I_b = params.ball_inertia;
R_eff = R_o - R_b; r = R_o/R_b;
a_old = m_b*R_eff^2 + I_b*r^2;
e_old = I_b*r*(r+1);
B4_before = e_old/a_old;
B4_after  = B_ref(4);
fprintf('  B(4) (u -> psi_ddot coupling):  BEFORE=%.4f  AFTER=%.4f  (%.1f%% change)  [%s]\n', ...
    B4_before, B4_after, 100*(B4_after-B4_before)/B4_before, pass_str(abs(B4_after-B4_before) > 1e-6));
[n_pass, n_fail] = tally(abs(B4_after-B4_before) > 1e-6, n_pass, n_fail);

% (b) plan_trajectory_casadi.m: symbolic dynamics now come from
% cascade_dynamics_reduced.m. Check with a real CasADi MX evaluation
% that xdot(4) (psi_ddot) matches the reference numerically at a test
% point, and differs from the old hand-derived formula.
import casadi.*
x_sym_t = MX.sym('x', 4);
u_sym_t = MX.sym('u', 1);
xdot_new = cascade_dynamics_reduced(x_sym_t, u_sym_t, params, 'rolling_out');
xdot_fun = Function('xdot_fun', {x_sym_t, u_sym_t}, {xdot_new});
x_num = [0; 0.3; deg2rad(5); 0.1];
u_num = 2.0;
psi_ddot_after = full(xdot_fun(x_num, u_num));
psi_ddot_after = psi_ddot_after(4);

% Old hand-derived formula from plan_trajectory_casadi.m (pre-fix), evaluated
% at the same point, using its own a/b/c/e (R_o-based, "+1" coupling):
b_old = params.ball_friction * r^2;
c_old = m_b * params.gravity * R_eff * sin(x_num(3));
psi_ddot_before = -(b_old*(x_num(4)-x_num(2)) + c_old + e_old*u_num) / a_old;

ok = abs(psi_ddot_after - psi_ddot_before) > 1e-6;
fprintf('  psi_ddot at test point:  BEFORE(old hand-derived)=%.4f  AFTER(cascade_dynamics_reduced)=%.4f  [%s]\n', ...
    psi_ddot_before, psi_ddot_after, pass_str(ok));
[n_pass, n_fail] = tally(ok, n_pass, n_fail);

fprintf('\n');

%% ---------------------------------------------------------------------
%  A4/Phase2 -- ball_hoop_ode.m returns a sol struct with mode history
%  as (t_start, t_end, mode) intervals
%% ---------------------------------------------------------------------
fprintf('--- A4/Phase2: ball_hoop_ode.m sol.mode_intervals ---\n');

% 'Free fall through the hole' (scenarios/open_loop/OL3_fall_through_gap.m):
% starts free_fall, should stay free_fall (falls through the gap, no
% rolling contact reached within the horizon for this IC) -- a simple
% case to check the interval list covers the run and starts in the
% right mode. Values inlined rather than calling OL3 directly: that
% file is a declaration that calls run_scenario(scn) (full integration
% + save + plot), not a bare (x0, mode) pair.
x0_3 = [params.outer_hoop_inner_radius - params.ball_radius - 0.001; 0; 0; 0; ...
        params.hole_center_angle + 0.15; 0];
state_3 = 'free_fall';
tau_fun_zero = @(t, x) 0;
sol = ball_hoop_ode( ...
    x0_3, tau_fun_zero, @ball_hoop_dynamics, params, cfg, state_3);

ok1 = isfield(sol, 't') && isfield(sol, 'X') && isfield(sol, 'te') && ...
      isfield(sol, 'Xe') && isfield(sol, 'mode_intervals');
fprintf('  sol has fields t/X/te/Xe/mode_intervals:  [%s]\n', pass_str(ok1));
[n_pass, n_fail] = tally(ok1, n_pass, n_fail);

ok2 = strcmp(sol.mode_intervals(1).mode, state_3);
fprintf('  mode_intervals(1).mode == initial_mode (''%s''):  [%s]\n', state_3, pass_str(ok2));
[n_pass, n_fail] = tally(ok2, n_pass, n_fail);

ok3 = abs(sol.mode_intervals(1).t_start - sol.t(1)) < 1e-15 && ...
      abs(sol.mode_intervals(end).t_end - sol.t(end)) < 1e-15;
fprintf('  intervals span exactly [t(1), t(end)]:  [%s]\n', pass_str(ok3));
[n_pass, n_fail] = tally(ok3, n_pass, n_fail);

% Every sample time falls within exactly one interval, and consecutive
% intervals are contiguous (no gaps, no overlaps).
ok4 = true;
for i = 2:numel(sol.mode_intervals)
    if abs(sol.mode_intervals(i).t_start - sol.mode_intervals(i-1).t_end) > 1e-12
        ok4 = false;
    end
end
fprintf('  intervals are contiguous (no gaps/overlaps):  [%s]\n', pass_str(ok4));
[n_pass, n_fail] = tally(ok4, n_pass, n_fail);

fprintf('\n');

%% ---------------------------------------------------------------------
%  A5 -- detect_events.m degenerate-state guard matches each mode's
%  declared event count
%% ---------------------------------------------------------------------
fprintf('--- A5: detect_events.m degenerate-guard event count ---\n');

x_degenerate = [0; 0; 0; 0; 0; 0];   % ball_r=0 < params.ball_radius, triggers the guard

modes_1event = {'rolling_out', 'rolling_in_outside'};
modes_2event = {'free_fall', 'rolling_in_inside'};

all_ok = true;
for k = 1:numel(modes_1event)
    [value, isterminal, direction] = detect_events(0, x_degenerate, params, cfg, modes_1event{k});
    ok = numel(value) == 1 && numel(isterminal) == 1 && numel(direction) == 1;
    fprintf('  %-20s degenerate guard: numel(value)=%d (expected 1, was 2 before)  [%s]\n', ...
        modes_1event{k}, numel(value), pass_str(ok));
    all_ok = all_ok && ok;
end
for k = 1:numel(modes_2event)
    [value, isterminal, direction] = detect_events(0, x_degenerate, params, cfg, modes_2event{k});
    ok = numel(value) == 2 && numel(isterminal) == 2 && numel(direction) == 2;
    fprintf('  %-20s degenerate guard: numel(value)=%d (expected 2, unaffected)  [%s]\n', ...
        modes_2event{k}, numel(value), pass_str(ok));
    all_ok = all_ok && ok;
end
[n_pass, n_fail] = tally(all_ok, n_pass, n_fail);

fprintf('\n');

%% ---------------------------------------------------------------------
%  Style 1 -- motor_inner_loop.m error() string-concat bug
%% ---------------------------------------------------------------------
fprintf('--- Style 1: motor_inner_loop.m error message construction ---\n');

% 'first_order' is now implemented (docs/MODEL.md sec. 10.5) and
% requires the tau_lag state as a 5th argument; calling without it is
% still an error, but a different, intentional one -- this check now
% confirms THAT error message is clean (not a string-concat/dimensions
% bug), not that the model itself is unimplemented.
x_phys_t = [0; 0; 0; 1; 0; 0];
cfg_first_order = cfg; cfg_first_order.motor_model = 'first_order';
try
    motor_inner_loop(1.0, x_phys_t, params, cfg_first_order);
    ok1 = false;
    msg1 = '(no error thrown)';
catch ME
    msg1 = ME.message;
    ok1 = contains(msg1, 'tau_lag') && ~contains(msg1, 'dimensions');
end
fprintf('  ''first_order'' (no tau_lag) error message: "%s"\n  [%s]\n', msg1, pass_str(ok1));
[n_pass, n_fail] = tally(ok1, n_pass, n_fail);

try
    tau3 = motor_inner_loop(1.0, x_phys_t, params, cfg_first_order, 0.1);
    ok1b = isfinite(tau3);
catch ME
    ok1b = false;
    fprintf('  ''first_order'' (with tau_lag) threw: %s\n', ME.message);
end
fprintf('  ''first_order'' (with tau_lag) runs without error:  [%s]\n', pass_str(ok1b));
[n_pass, n_fail] = tally(ok1b, n_pass, n_fail);

cfg_bogus = cfg; cfg_bogus.motor_model = 'bogus';
try
    motor_inner_loop(1.0, x_phys_t, params, cfg_bogus);
    ok2 = false;
    msg2 = '(no error thrown)';
catch ME
    msg2 = ME.message;
    ok2 = contains(msg2, 'bogus') && ~contains(msg2, 'dimensions');
end
fprintf('  unknown-model error message:  "%s"\n  [%s]\n', msg2, pass_str(ok2));
[n_pass, n_fail] = tally(ok2, n_pass, n_fail);

fprintf('\n');

%% ---------------------------------------------------------------------
%  Style 2 -- build_run.m/save_run.m use datetime instead of now/datestr
%% ---------------------------------------------------------------------
fprintf('--- Style 2: datetime instead of now/datestr ---\n');

scn_t = struct('name', 'validate_style2_test', ...
    'controller', struct('type', 'LQR'), ...
    'reference',  reference_constant(0));
t_t   = (0:0.1:1)';
X_t   = zeros(numel(t_t), 6);
tau_t = zeros(numel(t_t), 1);
u_t   = zeros(numel(t_t), 1);
run_t = build_run(t_t, X_t, tau_t, u_t, scn_t, params, cfg);

ok1 = isa(run_t.timestamp, 'datetime');
fprintf('  build_run(...).timestamp is a datetime object:  [%s]\n', pass_str(ok1));
[n_pass, n_fail] = tally(ok1, n_pass, n_fail);

tmp_folder = tempname;
filepath_t = save_run(run_t, tmp_folder);
[~, fname_t, fext_t] = fileparts(filepath_t);
fname_t = [fname_t fext_t];
ok2 = ~isempty(regexp(fname_t, '^\d{2}_\d{2}_\d{4}-\d{2}-\d{2}-\d{2}_lqr_validate_style2_test\.mat$', 'once'));
fprintf('  save_run(...) filename matches expected pattern: "%s"  [%s]\n', fname_t, pass_str(ok2));
[n_pass, n_fail] = tally(ok2, n_pass, n_fail);
if exist(tmp_folder, 'dir'), rmdir(tmp_folder, 's'); end

fprintf('\n');

%% ---------------------------------------------------------------------
%  Phase 2 -- run struct schema (results/RUN_SCHEMA.md), all fields
%  present with units in the name, and plot_closed_loop_results.m/
%  plot_runs.m read the real schema instead of a struct shape
%  build_run.m never produced (docs/AUDIT.md item C5)
%% ---------------------------------------------------------------------
fprintf('--- Phase 2: run struct schema ---\n');

% Minimal inline LQR scenario (shape: docs/MODEL.md sec. 9), not
% scenarios/closed_loop/T2/T2_tracking.m: that file always uses the full cfg (its
% point is to be runnable standalone), and this check wants a short
% horizon purely for speed.
cfg_short2 = cfg; cfg_short2.t_end = 2.0;
[~, R_eff_o] = hoop_geometry('rolling_out', params);
scn_schema_test.id   = 'validate_physics_schema_test';
scn_schema_test.name = 'validate_physics schema test';
scn_schema_test.mode = 'rolling_out';
scn_schema_test.x0   = [R_eff_o; 0; 0; 0; 0; 0];
scn_schema_test.controller = struct('type', 'LQR', 'psi_lin', 0, ...
    'Q', diag([1e-8, 1e-6, 1/deg2rad(3)^2, 1/0.3^2]), 'R', 5e-5);
scn_schema_test.reference  = reference_sinusoid(deg2rad(10), 0.15);

before_files = dir('runs/*.mat');
run_schema_test = run_scenario(scn_schema_test, params, cfg_short2);
close all;
after_files = dir('runs/*.mat');
new_names = setdiff({after_files.name}, {before_files.name});

expected_fields = {'name','t_s','x','tau_Nm','u_rad_s2','ref_psi_rad', ...
    'ref_psi_dot_rad_s','error_rad','z_ctrl','scenario','params','cfg','timestamp'};
ok = all(isfield(run_schema_test, expected_fields));
missing = expected_fields(~isfield(run_schema_test, expected_fields));
fprintf('  run struct has all schema fields:  [%s]', pass_str(ok));
if ~ok, fprintf('  (missing: %s)', strjoin(missing, ', ')); end
fprintf('\n');
[n_pass, n_fail] = tally(ok, n_pass, n_fail);

ok2 = size(run_schema_test.x, 2) == 6 && ...
      numel(run_schema_test.t_s) == size(run_schema_test.x, 1) && ...
      numel(run_schema_test.tau_Nm) == numel(run_schema_test.t_s) && ...
      numel(run_schema_test.u_rad_s2) == numel(run_schema_test.t_s);
fprintf('  array fields are consistently sized:  [%s]\n', pass_str(ok2));
[n_pass, n_fail] = tally(ok2, n_pass, n_fail);

% u_rad_s2 is genuinely new data (not previously logged): sanity-check
% it's not trivially all-zero for a tracking scenario that requires
% nonzero control effort.
ok3 = any(run_schema_test.u_rad_s2 ~= 0);
fprintf('  u_rad_s2 contains nonzero commanded accelerations:  [%s]\n', pass_str(ok3));
[n_pass, n_fail] = tally(ok3, n_pass, n_fail);

% plot_closed_loop_results.m and plot_runs.m read the schema without error.
try
    plot_closed_loop_results(run_schema_test, params);
    close(gcf);
    ok4 = true;
catch ME
    ok4 = false;
    fprintf('  plot_closed_loop_results threw: %s\n', ME.message);
end
fprintf('  plot_closed_loop_results(run, params) runs without error:  [%s]\n', pass_str(ok4));
[n_pass, n_fail] = tally(ok4, n_pass, n_fail);

try
    plot_runs(struct('files', {new_names(1)}));
    close all;
    ok5 = true;
catch ME
    ok5 = false;
    fprintf('  plot_runs threw: %s\n', ME.message);
end
fprintf('  plot_runs(...) runs without error (was broken pre-Phase-2, item C5):  [%s]\n', pass_str(ok5));
[n_pass, n_fail] = tally(ok5, n_pass, n_fail);

% plot_run.m (Phase 2 target item 6): single entry point, dispatches to
% plot_closed_loop_results.m / plot_states.m / plot_phase_portrait.m.
try
    plot_run(run_schema_test);                                            % default opts
    close all;
    plot_run(run_schema_test, struct('show_summary', false, ...
        'show_states', true, 'show_phase_portrait', true));               % all panels
    close all;
    ok6 = true;
catch ME
    ok6 = false;
    fprintf('  plot_run threw: %s\n', ME.message);
end
fprintf('  plot_run(run, opts) runs without error (default and all-panels opts):  [%s]\n', pass_str(ok6));
[n_pass, n_fail] = tally(ok6, n_pass, n_fail);

delete(fullfile('results', 'runs', new_names{1}));

fprintf('\n');

%% ---------------------------------------------------------------------
%  Summary
%% ---------------------------------------------------------------------
fprintf('==================================================================\n');
fprintf('TOTAL: %d passed, %d failed\n', n_pass, n_fail);
fprintf('==================================================================\n');

if n_fail > 0
    error('validate_physics: %d check(s) failed.', n_fail);
end


%% =======================================================================
%  Helpers
%% =======================================================================
function s = pass_str(ok)
%PASS_STR  Formats a boolean as 'PASS' or 'FAIL'.
    if ok, s = 'PASS'; else, s = 'FAIL'; end
end

function [np, nf] = tally(ok, np, nf)
%TALLY  Increments the pass or fail counter based on ok.
    if ok, np = np + 1; else, nf = nf + 1; end
end
