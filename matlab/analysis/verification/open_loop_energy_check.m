function results = open_loop_energy_check()
%OPEN_LOOP_ENERGY_CHECK  Mechanical energy along the open-loop scenarios,
%                        phase by phase and across every transition.
%
%   results = open_loop_energy_check()
%
%   Usage: matlab -batch "open_loop_energy_check"
%
%   WHAT IT TESTS, AND WHY IT IS NOT analysis/verification/run_all_checks.m CHECK 1.
%   Check 1 there conserves energy inside ONE rolling mode, which is a
%   test of the equations of motion and of the integrator. The claim the
%   thesis makes about the hybrid implementation is larger, and has three
%   parts, only the first of which check 1 covers:
%
%     (i)   within a continuous phase, E is constant;
%     (ii)  across a LIFTOFF, E is constant too -- no impulsive force acts
%           when contact is merely lost, so a jump there would mean the
%           reset map is inventing or destroying energy;
%     (iii) across an IMPACT, E may only DECREASE -- the landing map is
%           inelastic and imposes rolling instantly, and both of those
%           dissipate.
%
%   Parts (ii) and (iii) can only be seen on a run that actually leaves
%   and regains a surface, which is what the scenarios below add.
%
%   THE FREE-FLIGHT SPIN IS PART OF THE ENERGY. While the ball rolls, its
%   spin is fixed by the no-slip constraint and is already inside the
%   rolling M matrix. In flight it is a frozen constant that no state
%   carries (ball_spin_store.m), so it has to be reinstated here from the
%   last rolling state before liftoff -- otherwise E appears to jump at
%   every liftoff purely because a term was dropped from the bookkeeping.
%   This is also why analysis/verification/compute_energy.m is not
%   used: it accounts for spin on the outer track only.
%
%   Friction is switched off (ball_friction = motor_friction = 0) and no
%   torque is applied, so the only mechanisms that may change E are the
%   impact map and integration error.
%
%   Outputs
%   -------
%   results : struct array, one entry per scenario, with fields
%               .scenario     scenario id
%               .phases       struct array (.mode, .t_start, .t_end,
%                             .drift_pct) -- worst deviation from the
%                             phase's own initial energy, in % of |E(0)|
%               .transitions  struct array (.t, .from, .to, .jump_J,
%                             .jump_pct_T) -- the energy step across the
%                             reset map, in joules and as a percentage of
%                             the KINETIC energy just before it. The
%                             kinetic energy is the reference here, not
%                             E(0): the potential is measured from the
%                             hoop axis, so E(0) can be small or of either
%                             sign, and a step reported against it says
%                             more about where the datum sits than about
%                             the impact.
%               .E0           initial total mechanical energy [J]

    project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(genpath(project_root));

    % The scenarios are named, never restated: their initial conditions
    % live in scenarios/open_loop/*.m and are read back from the run each
    % one produces, so this check cannot drift out of step with them.
    scenario_names = { ...
        'OL1_rolling_outer', ...
        'OL2_rolling_inner_inside', ...
        'OL9_escape_inner_hoop', ...
        'OL10_energy_validation'};

    T_END = 3.0;   % [s] horizon the thesis quotes for the energy figures

    params = ball_hoop_params();
    params.ball_friction  = 0;
    params.motor_friction = 0;

    cfg = sim_config();
    cfg.t_end = T_END;

    fprintf('\n==========================================================\n');
    fprintf('OPEN_LOOP_ENERGY_CHECK  (tau = 0, friction off, t_end = %.1f s)\n', T_END);
    fprintf('==========================================================\n');

    results = struct('scenario', {}, 'phases', {}, 'transitions', {}, 'E0', {});

    for k = 1:numel(scenario_names)
        name = scenario_names{k};

        % Run the scenario as declared, only to recover its initial
        % condition and mode; the run itself (default params, with
        % friction) is discarded.
        declared = feval(name);
        close all;
        x0   = declared.scenario.x0;
        mode = declared.scenario.mode;

        ball_spin_store('reset');
        sol = ball_hoop_ode(x0, @(t, x) 0, @ball_hoop_dynamics, params, cfg, mode);

        r = energy_report(sol, params);
        r.scenario = name;
        results(end+1) = r; %#ok<AGROW>

        print_report(r);
    end

    close all;
end


function r = energy_report(sol, params)
%ENERGY_REPORT  Splits a solution into its continuous phases, evaluates the
%               total mechanical energy on each, and measures the drift
%               inside each phase and the jump across each transition.
    t  = sol.t;
    X  = sol.X;
    mi = sol.mode_intervals;

    E   = zeros(numel(t), 1);
    T   = zeros(numel(t), 1);
    phases = struct('mode', {}, 't_start', {}, 't_end', {}, 'drift_pct', {});

    % Spin the ball carries into a flight: set at the last rolling sample
    % before liftoff, exactly as transition_state.m does, and zero if the
    % scenario starts airborne (a ball released from rest is not spinning).
    phi_dot_flight = 0;
    idx_last = [];   % index of the final sample of the previous phase

    for j = 1:numel(mi)
        if j == 1
            idx = find(t >= mi(j).t_start & t <= mi(j).t_end);
        else
            idx = find(t > mi(j).t_start & t <= mi(j).t_end);
        end
        if isempty(idx), continue; end

        for i = idx(:)'
            [E(i), T(i)] = total_energy(X(i, :).', mi(j).mode, phi_dot_flight, params);
        end

        phases(end+1) = struct('mode', mi(j).mode, ...
            't_start', mi(j).t_start, 't_end', mi(j).t_end, ...
            'drift_pct', max(abs(E(idx) - E(idx(1)))) / abs(E(idx(1))) * 100); %#ok<AGROW>

        if ~strcmp(mi(j).mode, 'free_fall')
            % Valid until the ball leaves this surface, which is exactly
            % when the value is needed.
            phi_dot_flight = ball_spin_from_rolling(X(idx(end), :).', params, mi(j).mode);
        end

        idx_last(end+1) = idx(end); %#ok<AGROW>
    end

    transitions = struct('t', {}, 'from', {}, 'to', {}, 'jump_J', {}, 'jump_pct_T', {});
    for j = 2:numel(phases)
        i_pre  = idx_last(j-1);
        i_post = find(t > mi(j).t_start, 1, 'first');
        transitions(end+1) = struct('t', phases(j).t_start, ...
            'from', phases(j-1).mode, 'to', phases(j).mode, ...
            'jump_J', E(i_post) - E(i_pre), ...
            'jump_pct_T', (E(i_post) - E(i_pre)) / T(i_pre) * 100); %#ok<AGROW>
    end

    r = struct('scenario', '', 'phases', phases, 'transitions', transitions, 'E0', E(1));
end


function [E, T] = total_energy(x, mode, phi_dot_flight, params)
%TOTAL_ENERGY  E = T + V for one state, in the mode that state belongs to,
%              returning the kinetic part separately. Potential is
%              referenced to the hoop axis, so it is negative at the bottom
%              of the track (psi = 0).
    m = params.ball_mass;
    g = params.gravity;

    if strcmp(mode, 'free_fall')
        r_ball  = x(1);
        r_dot   = x(2);
        psi     = x(5);
        psi_dot = x(6);
        T = 0.5 * m * (r_dot^2 + (r_ball * psi_dot)^2) ...
          + 0.5 * params.hoop_motor_inertia * x(4)^2 ...
          + 0.5 * params.ball_inertia * phi_dot_flight^2;
        V = - m * g * r_ball * cos(psi);
    else
        [~, R_eff] = hoop_geometry(mode, params);
        [M, ~, ~]  = rolling_matrices(params, x(5), mode);
        qdot = [x(4); x(6)];
        T = 0.5 * qdot' * M * qdot;
        V = - m * g * R_eff * cos(x(5));
    end
    E = T + V;
end


function print_report(r)
    fprintf('\n%s   (E(0) = %+.6f J)\n', r.scenario, r.E0);
    fprintf('  phases:\n');
    for j = 1:numel(r.phases)
        p = r.phases(j);
        fprintf('    [%7.4f %7.4f] s  %-18s  max drift %.4f %%\n', ...
            p.t_start, p.t_end, p.mode, p.drift_pct);
    end
    if isempty(r.transitions)
        fprintf('  transitions: none\n');
        return;
    end
    fprintf('  transitions:\n');
    for j = 1:numel(r.transitions)
        tr = r.transitions(j);
        if strcmp(tr.from, 'free_fall')
            kind = 'impact ';
        else
            kind = 'liftoff';
        end
        fprintf('    t = %7.4f s  %s  %-18s -> %-18s  dE = %+.3e J  (%+.2f %% of T-)\n', ...
            tr.t, kind, tr.from, tr.to, tr.jump_J, tr.jump_pct_T);
    end
end
