function [t_traj, x_traj, u_traj, Tf_opt, sol] = ...
    plan_trajectory_casadi(x0, xf, opts, params, mode)
%PLAN_TRAJECTORY_CASADI  CasADi/IPOPT direct-collocation trajectory planner.
%
%   Inputs
%   ------
%   x0, xf : boundary reduced states [theta; theta_dot; psi; psi_dot]
%            [rad; rad/s; rad; rad/s]
%   opts   : struct with fields:
%              .K      : number of collocation points
%              .Tf_min : minimum trajectory duration [s]
%              .Tf_max : maximum trajectory duration [s]
%              .u_max  : control saturation limit [rad/s^2]
%              .contact_margin : (optional, default 0) lower bound on
%                        N/m = g*cos(psi) + R_eff*psi_dot^2, the normal
%                        force per unit ball mass [m/s^2]. Zero means
%                        "never pull on the ball", the bare physical
%                        no-detachment condition -- but since the cost
%                        minimises control effort and extra speed costs
%                        effort, the optimum then rides exactly on
%                        N = 0 and is untrackable in closed loop (any
%                        tracking error detaches the ball; verified for
%                        the loop-the-loop case across five decades of
%                        TVLQR R). A positive margin buys headroom: at
%                        the top of a loop it raises the required speed
%                        from sqrt(g/R_eff) to sqrt((g+margin)/R_eff),
%                        so margin = 0.1*g costs only ~5% more speed
%                        there. Expressing it as a fraction of
%                        params.gravity keeps it dimensionally readable.
%              .psi_bounds : (optional, default [-3.5*pi, 0.5*pi]) search
%                        box on psi [rad]. The default is sized for a
%                        revolution in the -psi direction; a manoeuvre
%                        that drives psi positive past +pi/2 must widen
%                        it, or the NLP is infeasible for a reason that
%                        has nothing to do with the physics.
%              .liftoff : (optional, default absent) struct with field
%                        .psi_range = [psi_min, psi_max] in rad. When
%                        present the plan ends AT DETACHMENT instead of
%                        at xf: the terminal knot is constrained to
%                        N = 0 with psi(Tf) free inside psi_range and
%                        psi_dot(Tf) >= 0, and xf is used only to seed
%                        the initial guess. psi_range must lie in
%                        (pi/2, pi] -- detachment needs cos(psi) < 0 --
%                        and inside psi_bounds. This is T4 phase 1.
%                        Optional extra fields:
%                          .theta_range     [min max] window on the RAW
%                            theta(Tf) [rad], i.e. on the hoop
%                            orientation at release. Windings are not
%                            wrapped, so the caller places the window at
%                            the intended multiple of 2*pi. Use it to
%                            guarantee the phase-2 repositioning is
%                            reachable instead of checking afterwards.
%                          .reach_radius    the ballistic path must pass
%                            within this radius of the hoop axis [m].
%                            Without it the optimum releases at the
%                            cheapest angle in the window and the ball
%                            goes nowhere -- see the constraint block.
%                          .max_flight_time upper bound on the free
%                            flight time used by that constraint [s],
%                            default 1.0.
%              .n_swings : (optional, default 0) number of pump-up
%                        swings written into the initial guess. 0 keeps
%                        the monotone cosine ramp on psi, which is the
%                        right hint for a one-shot manoeuvre and the
%                        wrong topological class for a loop-the-loop.
%                        See build_initial_guess.
%              .constrain_theta_f : (optional, default false) also pin
%                        theta(Tf) to xf(1). Off by default because a
%                        loop-the-loop does not care where the hoop ends
%                        up; T4 does, because the hoop carries the gap
%                        the ball has to fall through.
%              .warm_start : (optional, default absent) struct with
%                        fields .X (nx x K), .U (1 x K), .Tf (scalar) used
%                        as the initial guess in place of the
%                        resonance-pumping heuristic. This is the same
%                        mechanism the continuation fallback already uses
%                        internally between its rungs, exposed to callers
%                        that walk a parameter themselves -- a margin
%                        sweep re-solves the SAME boundary problem once
%                        per margin, and the previous converged plan is a
%                        far better starting point than the generic guess.
%                        The natural chaining is
%                            opts.warm_start = struct('X', x_traj.', ...
%                                                     'U', u_traj.', ...
%                                                     'Tf', Tf_opt);
%                        from the previous solve's outputs. Absent (or
%                        empty) reproduces the previous behaviour exactly,
%                        so no existing caller changes. Note that a warm
%                        start selects a BRANCH as much as it accelerates
%                        convergence -- see the continuation fallback's
%                        note on the 2.4 s and 4.5 s branches -- so a
%                        sweep warm-started along a parameter is tracking
%                        one branch continuously, which is usually what is
%                        wanted but is a different object from a set of
%                        independent cold solves.
%   params : physical parameters from ball_hoop_params()
%   mode   : 'rolling_out' | 'rolling_in_outside' | 'rolling_in_inside'
%            (canonical mode strings, see hoop_geometry.m -- this function
%            used to accept its own 'outer'/'inner_out'/'inner_in' names
%            and translate them internally; Phase 2 removed that
%            translation layer so there is exactly one mode vocabulary in
%            the codebase, not two. Callers updated accordingly.)
%
%   Outputs
%   -------
%   t_traj : time grid (K x 1)                              [s]
%   x_traj : optimized reduced state trajectory (K x 4)
%            [theta, theta_dot, psi, psi_dot], [rad, rad/s, rad, rad/s]
%   u_traj : optimized control trajectory (K x 1)             [rad/s^2]
%   Tf_opt : optimized trajectory duration                     [s]
%   sol    : CasADi Opti solution object (diagnostics, solver stats)
%
%   Robustness note: when opts.u_max is close to the true feasibility
%   boundary for (x0, xf) -- e.g. a fast full loop near the actuator's
%   real torque limit -- the single-shot resonance-pumping initial guess
%   below can fail to converge even though a feasible trajectory exists
%   (IPOPT reports "locally infeasible" from a poor starting point, not a
%   genuine infeasibility proof; verified by hand for the T3
%   loop-the-loop case). If the direct solve fails, this function
%   automatically retries via continuation on u_max: solve first at a
%   comfortably large torque budget, then bisect down toward the
%   requested opts.u_max, warm-starting each solve from the previous one.
%   Callers whose direct solve already converges see no behaviour change.

    % A caller-supplied warm start replaces the heuristic guess for the
    % DIRECT solve only. If that still fails, the continuation fallback
    % below starts over from its own easy point, exactly as before: a
    % warm start that did not work is not a better place to restart from.
    if isfield(opts, 'warm_start') && ~isempty(opts.warm_start)
        warm = opts.warm_start;
    else
        warm = [];
    end

    try
        [t_traj, x_traj, u_traj, Tf_opt, sol] = solve_once(x0, xf, opts, params, mode, warm);
    catch ME_direct
        fprintf(['\nDirect solve failed (u_max=%.3f); retrying via u_max continuation ' ...
                  'from a larger torque budget...\n'], opts.u_max);
        [t_traj, x_traj, u_traj, Tf_opt, sol] = solve_via_continuation(x0, xf, opts, params, mode, ME_direct);
    end
end


%% ========================================================================
%  CORE SOLVE (single u_max)
%% ========================================================================

function [t_traj, x_traj, u_traj, Tf_opt, sol, X_val, U_val, Tf_val] = ...
    solve_once(x0, xf, opts, params, mode, warm, report_on_failure)
%SOLVE_ONCE  Builds and solves the direct-collocation NLP for a single
%            opts.u_max value. warm, if non-empty, is a struct with
%            fields .X (nx x K), .U (1 x K), .Tf (scalar) used as the
%            initial guess instead of the resonance-pumping heuristic --
%            X_val/U_val/Tf_val (raw Opti values) are returned alongside
%            the usual outputs specifically so a caller can chain this as
%            the next warm start. report_on_failure (default true)
%            suppresses the diagnostic figure (console output stays) when
%            false -- used by the continuation fallback so an internal
%            bisection step that fails, as expected, doesn't pop a figure
%            for every attempt.

if nargin < 7 || isempty(report_on_failure)
    report_on_failure = true;
end

import casadi.*

%% Symbolic Dynamics
nx = 4;
nu = 1;

x_sym = MX.sym('x', nx);
u_sym = MX.sym('u', nu);

% Dynamics come from cascade_dynamics_reduced.m (itself built on
% rolling_matrices.m, the single source of truth for this system) --
% directly CasADi-compatible: it only uses +,-,*,/,sin on its inputs,
% all of which CasADi's MX class overloads. An earlier, independent
% hand-derivation here used R_o (hoop radius) where the ball's own-spin
% term needs R_eff (ball-centre orbit radius), and a "+1" where the
% coupling term needs "-1" -- both silently wrong (docs/AUDIT.md items
% 1.7, C6).
[~, R_eff, ~, ~, ~, r_roll] = hoop_geometry(mode, params);   % still needed below (contact constraint, initial guess)
g = params.gravity;                          % still needed below (contact constraint, initial guess)

% Two-rail rolling factor: the ball's effective inertia about the hoop
% axis is m*R_eff^2*k (rolling_matrices.m line 50 factorised), so the
% pendulum frequency of the ball on a HELD hoop is sqrt(g/(k*R_eff)).
% Only build_initial_guess uses it, but it is computed here because
% r_roll is already unpacked above.
k_roll = 1 + params.ball_inertia / (params.ball_mass * r_roll^2);

xdot = cascade_dynamics_reduced(x_sym, u_sym, params, mode);

f = Function('f', {x_sym, u_sym}, {xdot}, {'x','u'}, {'xdot'});             % xdot = f(x_sym, u_sym)

%% Optimization problem
K      = opts.K;
Tf_min = opts.Tf_min;
Tf_max = opts.Tf_max;
u_max  = opts.u_max;
if isfield(opts, 'contact_margin') && ~isempty(opts.contact_margin)
    contact_margin = opts.contact_margin;
else
    contact_margin = 0;   % bare no-detachment condition (see docstring)
end
constrain_theta_f = isfield(opts, 'constrain_theta_f') && ~isempty(opts.constrain_theta_f) ...
                    && opts.constrain_theta_f;
if isfield(opts, 'psi_bounds') && ~isempty(opts.psi_bounds)
    psi_bounds = opts.psi_bounds;
else
    % Default covers a full revolution in the -psi direction (the
    % loop-the-loop case) with margin beyond [-pi, 0]. A manoeuvre that
    % drives psi POSITIVE past +pi/2 -- T4's liftoff at +125 deg -- must
    % widen this or the NLP is infeasible for a reason that has nothing
    % to do with the physics.
    psi_bounds = [-3.5*pi, 0.5*pi];
end

% Velocity boxes. These were literals of 6*pi described as "arbitrary
% generous bounds, well beyond any trajectory this planner is expected to
% produce". That was false for the looping case: measured on T3b, the
% solution sits exactly ON the theta_dot bound (18.8496 = 6*pi to four
% decimals), so the plan was being shaped by a number chosen to be
% meaningless. They are options now, still defaulting to 6*pi so no
% existing scenario changes, and a caller with a real limit -- the
% ODrive's vel_limit, 20 rad/s (t1_config_v3.py, THETA_DOT_MAX) -- should
% pass it instead of inheriting this one.
if isfield(opts, 'theta_dot_max') && ~isempty(opts.theta_dot_max)
    theta_dot_max = opts.theta_dot_max;
else
    theta_dot_max = 6*pi;
end
if isfield(opts, 'psi_dot_max') && ~isempty(opts.psi_dot_max)
    psi_dot_max = opts.psi_dot_max;
else
    psi_dot_max = 6*pi;
end

% No-slip. OPT-IN, default off: turning it on changes every trajectory
% this planner has ever produced, and the two scenarios that predate it
% (T3b, T4b) are on record with the numbers they were measured at.
% Measured on T3b's plan, 17 of its 60 knots need more friction than the
% contact can supply -- and 10 of them still do at mu = 1, because the
% shortfall near the top is N falling toward zero, not mu being small.
% See dynamics/reset_maps/contact_forces.m.
% Number of pump-up swings written into the initial guess. 0 keeps the
% monotone cosine ramp every scenario before T3c was solved with, so
% those solves are unchanged; see build_initial_guess for what a
% positive value does and why the ramp is a problem for a looping plan.
if isfield(opts, 'n_swings') && ~isempty(opts.n_swings)
    n_swings = opts.n_swings;
else
    n_swings = 0;
end

% Liftoff terminal condition (T4 phase 1). When set, the plan does NOT
% end at a prescribed state: it ends wherever the ball detaches, and the
% detachment angle is a decision variable. See the constraint block below.
if isfield(opts, 'liftoff') && ~isempty(opts.liftoff)
    liftoff = opts.liftoff;
    if ~isfield(liftoff, 'psi_range') || numel(liftoff.psi_range) ~= 2
        error('plan_trajectory_casadi:liftoff_psi_range', ...
            'opts.liftoff needs a .psi_range = [psi_min, psi_max] in rad.');
    end
    if liftoff.psi_range(2) > psi_bounds(2) || liftoff.psi_range(1) < psi_bounds(1)
        error('plan_trajectory_casadi:liftoff_outside_box', ...
            ['opts.liftoff.psi_range = [%.3f, %.3f] is not inside opts.psi_bounds ' ...
             '= [%.3f, %.3f]; the terminal condition would be unreachable for a ' ...
             'reason that has nothing to do with the physics.'], ...
            liftoff.psi_range(1), liftoff.psi_range(2), psi_bounds(1), psi_bounds(2));
    end
else
    liftoff = [];
end

enforce_no_slip = isfield(opts, 'enforce_no_slip') && ~isempty(opts.enforce_no_slip) ...
                  && opts.enforce_no_slip;
constrain_midpoints = isfield(opts, 'constrain_midpoints') && ~isempty(opts.constrain_midpoints) ...
                      && opts.constrain_midpoints;

opti = Opti();

% Declare only K-1 free state nodes: column 1 is fixed to x0 by
% construction (see X = [x0, X_free] below), not as an optimization variable.
X_free = opti.variable(nx, K-1);
X = [x0, X_free];                   % column 1 is fixed by construction

U  = opti.variable(nu, K);
Tf = opti.variable();

h = Tf / (K-1);   % K-1 intervals: node 1 is fixed by construction, not a free variable

%% Collocation constraints (Hermite-Simpson)
% The midpoint states are kept so the physical constraints below can be
% imposed there too, not only at the knots: with K = 60 knots over a 2.4 s
% loop the gaps are 40 ms wide, long enough for a constraint satisfied at
% both ends to be violated in between.
Xc_cells = cell(1, K-1);
Uc_cells = cell(1, K-1);
for k = 1:K-1
    f_k   = f(X(:,k),   U(:,k));
    f_kp1 = f(X(:,k+1), U(:,k+1));

    x_c    = 0.5 * (X(:,k) + X(:,k+1)) + (h/8) * (f_k - f_kp1);
    xdot_c = -(1.5/h) * (X(:,k) - X(:,k+1)) - 0.25 * (f_k + f_kp1);

    u_c = 0.5 * (U(:,k) + U(:,k+1));

    opti.subject_to(f(x_c, u_c) - xdot_c == 0);

    Xc_cells{k} = x_c;
    Uc_cells{k} = u_c;
end
Xc = [Xc_cells{:}];
Uc = [Uc_cells{:}];

%% Boundary, control, physical constraints

% --- Final state: theta_dot, psi, psi_dot are always constrained ---
% theta is left free by default: for a loop-the-loop nothing depends on
% where the hoop ends up. T4 is the exception -- the hoop carries the
% inner-hoop gap, so the ball can only be caught if theta is in the right
% place at liftoff, and steering it there during the ~0.14 s of flight
% would need most of the torque budget. Constraining it here instead
% costs nothing and makes the flight phase a pure hold.
if isempty(liftoff)
    opti.subject_to(X(2:4, end) == xf(2:4));
else
    % T4 PHASE 1: the plan ends AT DETACHMENT, and where that happens is
    % for the optimiser to choose. Pinning it (as T4b does, at 132 deg)
    % forces a specific point out of a one-parameter family, on the
    % strength of a sweep run outside the optimisation. Letting psi_L
    % float inside a window instead lets the planner trade the release
    % point against the effort it costs to get there.
    %
    % Detachment IS N = 0, i.e. g*cos(psi) + R_eff*psi_dot^2 = 0 (see the
    % contact block below, where this is imposed on the terminal knot).
    % So the release SPEED is not an extra degree of freedom: it follows
    % from the angle, psi_dot_L = sqrt(-g*cos(psi_L)/R_eff), and only
    % exists at all for cos(psi_L) < 0 -- the upper half of the hoop,
    % which is why the window must live in (pi/2, pi].
    %
    % psi_dot(end) >= 0 picks the branch: the ball must be travelling UP
    % the hoop when it lets go. Without it the same N = 0 condition is
    % satisfied by a ball falling back down through the mirror state,
    % which is a detachment but not a launch.
    opti.subject_to(liftoff.psi_range(1) <= X(3,end) <= liftoff.psi_range(2));
    opti.subject_to(X(4,end) >= 0);

    % HOOP ORIENTATION AT RELEASE, as a window rather than a point.
    % Leaving theta free makes it an accident of whichever branch IPOPT
    % lands on -- MEASURED, the same 125 deg release ends at theta_L =
    % -106 deg or +125 deg depending only on how the psi window was
    % written -- and phase 2 then has to live with whatever it gets.
    % Pinning it to a single value (T4b's approach) over-constrains: it
    % forces the ball to coast the last stretch while the controller
    % holds the hoop, and holding the hoop against the rolling ball
    % pumps energy through the M(1,2) coupling.
    %
    % A window is the middle course, and it is the same device the
    % contact margin uses on N: state the range that keeps the next
    % phase feasible, and let the optimiser move inside it.
    %
    % NOTE ON WINDING. theta is not wrapped anywhere in this planner --
    % it accumulates, and a loop-the-loop plan ends at several turns. So
    % this bound is on the RAW theta and therefore fixes the number of
    % turns as well as the orientation. Only orientation matters
    % physically, so a caller wanting "orientation in [a, b], any
    % winding" must place the window at the right multiple of 2*pi
    % itself; t4c_phase1_liftoff.m does that from a first free solve.
    if isfield(liftoff, 'theta_range') && ~isempty(liftoff.theta_range)
        if numel(liftoff.theta_range) ~= 2 || diff(liftoff.theta_range) <= 0
            error('plan_trajectory_casadi:liftoff_theta_range', ...
                'opts.liftoff.theta_range must be [theta_min, theta_max] with theta_max > theta_min.');
        end
        opti.subject_to(liftoff.theta_range(1) <= X(1,end) <= liftoff.theta_range(2));
    end

    % theta_dot at liftoff is deliberately left free inside its box. It
    % is the phase-2 lever: once the ball is airborne the hoop decouples
    % completely (free_fall_matrices.m does not couple theta to r, psi),
    % so whatever speed the hoop carries into the flight is free
    % repositioning toward the catch angle.

    % REACHABILITY. Without this, minimising effort with a free release
    % point has an obvious and useless optimum: detach at the lowest
    % angle the window allows, where N = 0 costs the least speed, and go
    % nowhere. MEASURED on a [95, 175] deg window with no reach
    % constraint: the optimiser picks exactly 95.000 deg and the ball's
    % closest approach to the axis is 89.2 mm, against the 32.6 mm it
    % needs to get inside the inner hoop.
    %
    % The fix is to make phase 1 answer for what phase 2 needs, which is
    % the coupling Betts ch. 4 argues for rather than a sequential solve.
    % Once N = 0 the ball is ballistic and the hoop is irrelevant to WHERE
    % it goes, so the flight is closed form: with psi measured from the
    % downward vertical and the frame x-right / y-up, the ball centre
    % leaves from R_eff*[sin psi; -cos psi] with the tangential velocity
    % R_eff*psi_dot*[cos psi; sin psi].
    %
    % "The path passes within reach_radius of the axis" is imposed via a
    % free flight time rather than by minimising the distance over t:
    % EXISTS t >= 0 such that |p(t)| <= reach_radius is exactly the same
    % statement, and it is smooth, where the minimising t solves a cubic.
    if isfield(liftoff, 'reach_radius') && ~isempty(liftoff.reach_radius)
        if isfield(liftoff, 'max_flight_time') && ~isempty(liftoff.max_flight_time)
            t_fly_max = liftoff.max_flight_time;
        else
            t_fly_max = 1.0;   % far longer than any flight this geometry produces
        end
        t_fly = opti.variable();
        opti.subject_to(0 <= t_fly <= t_fly_max);

        psi_L     = X(3,end);
        psi_dot_L = X(4,end);
        px = R_eff*sin(psi_L)  + R_eff*psi_dot_L*cos(psi_L)*t_fly;
        py = -R_eff*cos(psi_L) + R_eff*psi_dot_L*sin(psi_L)*t_fly - 0.5*g*t_fly^2;

        opti.subject_to(px^2 + py^2 <= liftoff.reach_radius^2);
        opti.set_initial(t_fly, 0.15);
    end
end
if constrain_theta_f
    opti.subject_to(X(1, end) == xf(1));
end

% --- Control bounds ---
opti.subject_to(-u_max <= U <= u_max);

% --- Final time bounds ---
% Tf is itself a free decision variable, so it needs explicit bounds.
opti.subject_to(Tf_min <= Tf <= Tf_max);

% --- Physical bounds on states ---
% Theta unbounded
% +/-6*pi (3 revolutions) is an arbitrary generous bound, well beyond
% any trajectory this planner is expected to produce -- it exists only
% so the NLP has a finite search box, not as a physical limit.
opti.subject_to(-theta_dot_max <= X(2,:) <= theta_dot_max);
opti.subject_to(psi_bounds(1) <= X(3,:) <= psi_bounds(2));
opti.subject_to(-psi_dot_max <= X(4,:) <= psi_dot_max);

% --- Contact constraint ---
% Normal force N = m*(g*cos(psi) + R_eff*psi_dot^2) must stay non-negative
% (see plot_normal_force.m for the same free-body relation), or above a
% positive contact_margin when the caller asks for tracking headroom.
% docs/MODEL.md sec. 6.5 describes this as imposed "at every knot and
% midpoint"; the code only ever did the knots. Adding the midpoints is a
% real tightening -- measured, it moves T3b's optimum from Tf = 2.38 s to
% 3.23 s and breaks T4b's catch outright -- so it is opt-in rather than
% retrofitted onto two scenarios whose published numbers were measured
% without it. New work should set it; docs/MODEL.md has been corrected to
% say which scenarios do.
if isempty(liftoff)
    opti.subject_to(g * cos(X(3,:)) + R_eff * X(4,:).^2 >= contact_margin);
    if constrain_midpoints
        opti.subject_to(g * cos(Xc(3,:)) + R_eff * Xc(4,:).^2 >= contact_margin);
    end
else
    % A tracking margin and a liftoff terminal condition are asking for
    % opposite things, so they are separated in TIME rather than traded
    % off: the margin holds everywhere the ball is still meant to be
    % attached, and the LAST knot is the release, where N is exactly
    % zero. The two are not in conflict across one interval -- at the
    % release angle d(N/m)/dt is of order -g*sin(psi)*psi_dot ~ -48 m/s^3,
    % so a 40-70 ms interval crosses a 0.1*g margin comfortably.
    opti.subject_to(g * cos(X(3,1:end-1)) + R_eff * X(4,1:end-1).^2 >= contact_margin);
    opti.subject_to(g * cos(X(3,end))     + R_eff * X(4,end)^2      == 0);
    if constrain_midpoints
        opti.subject_to(g * cos(Xc(3,1:end-1)) + R_eff * Xc(4,1:end-1).^2 >= contact_margin);
        % The final midpoint sits between the last margin-carrying knot
        % and the release itself, so it may only be required to keep
        % contact, not to keep headroom.
        opti.subject_to(g * cos(Xc(3,end)) + R_eff * Xc(4,end)^2 >= 0);
    end
end

% --- No-slip constraint ---
% Staying ON the track (N >= 0 above) is not the same as ROLLING on it.
% The friction the contact must supply to keep the ball rolling is capped
% at mu*N, and near the top of a loop N is small, so the cap collapses
% exactly where a minimum-effort planner is most tempted to spend control.
% dynamics/reset_maps/contact_forces.m derives both forces and is checked against
% rolling_matrices.m row 2 numerically (analysis/verification/run_all_checks.m).
if enforce_no_slip
    [~, ~, slip_expr] = contact_forces(x_sym, u_sym, params, mode);
    slip_fun = Function('slip', {x_sym, u_sym}, {slip_expr});
    for k = 1:K
        opti.subject_to(slip_fun(X(:,k), U(:,k)) >= 0);
    end
    for k = 1:K-1
        opti.subject_to(slip_fun(Xc(:,k), Uc(:,k)) >= 0);
    end
end

%% Initial guess: warm-started from a previous solve, or resonance pumping
if isempty(warm)
    % A liftoff plan ends MOVING, at the speed N = 0 implies, so the guess
    % must end moving too: the zero-slope ramp otherwise hands IPOPT a
    % starting point that violates the terminal condition by the entire
    % release speed, and the solve fails from a bad guess rather than
    % from an infeasible problem (measured: it does).
    if isempty(liftoff)
        psi_dot_f_guess = [];
    else
        psi_dot_f_guess = xf(4);
    end
    [X_init, U_init, Tf_init] = build_initial_guess(f, x0, xf, K, Tf_min, Tf_max, R_eff, g, u_max, k_roll, n_swings, psi_dot_f_guess);
else
    X_init = warm.X;
    X_init(:,1) = x0;                          % column 1 must exactly match the fixed first column
    U_init = min(max(warm.U, -u_max), u_max);   % clip to this step's (possibly tighter) bound
    Tf_init = min(max(warm.Tf, Tf_min), Tf_max);
end
opti.set_initial(Tf, Tf_init);
opti.set_initial(U, U_init);
opti.set_initial(X, X_init);

%% Cost function
J = 0;
for k = 1:K-1
    U_k    = U(:,k);
    U_kp1  = U(:,k+1);
    U_c    = 0.5 * (U_k + U_kp1);
    J = J + (h/6) * (U_k.^2 + 4*U_c.^2 + U_kp1.^2);
end

opti.minimize(J);


%% IPOPT and solve
% Iteration budget and tolerances below are generous, arbitrary defaults
% (not derived from a convergence study for this problem size).
ipopt_opts = struct();
ipopt_opts.print_level             = 5;
ipopt_opts.max_iter                = 3000;
ipopt_opts.tol                     = 1e-6;
ipopt_opts.constr_viol_tol         = 1e-6;
ipopt_opts.acceptable_tol          = 1e-4;
ipopt_opts.hessian_approximation   = 'limited-memory';

opti.solver('ipopt', struct(), ipopt_opts);

try
    sol = opti.solve();
catch ME
    if report_on_failure
        report_failed_solve(opti, X, U, Tf, K, opts.u_max);
    end
    rethrow(ME);
end

%% Extract results
x_traj = sol.value(X)';
u_traj = sol.value(U)';
Tf_opt = sol.value(Tf);
t_traj = linspace(0, Tf_opt, K)';
X_val  = sol.value(X);
U_val  = sol.value(U);
Tf_val = Tf_opt;

fprintf('\nSolve succeeded.\n');
fprintf('  Tf_opt        = %.4f s\n', Tf_opt);
fprintf('  max |u|       = %.2f rad/s^2 (bound: %.2f)\n', max(abs(u_traj)), u_max);
fprintf('  theta_final   = %.2f rad (= %.2f turns)\n', x_traj(end,1), x_traj(end,1)/(2*pi));
fprintf('  psi_final     = %.2f deg (target: %.2f)\n', ...
        rad2deg(x_traj(end,3)), rad2deg(xf(3)));
fprintf('  cost J        = %.4e\n', sol.value(J));

end


%% ========================================================================
%  CONTINUATION FALLBACK
%% ========================================================================

function [t_traj, x_traj, u_traj, Tf_opt, sol] = ...
    solve_via_continuation(x0, xf, opts, params, mode, ME_direct)
%SOLVE_VIA_CONTINUATION  Fallback used by plan_trajectory_casadi.m when
%   the direct, single-shot solve fails. Solves first at a comfortably
%   large u_max (known to be far from any feasibility boundary), then
%   walks the torque budget down toward the originally-requested
%   opts.u_max through a GEOMETRIC ladder of intermediate steps (not a
%   single large jump), warm-starting each solve from the previous one.
%
%   Gradual steps matter here, not just eventual convergence: this
%   problem has multiple local solution branches (verified by hand for
%   the T3 loop-the-loop case -- a ~4.5s branch that reaches the real
%   torque budget, and a faster ~2.4s branch that does not). Jumping
%   straight from the easy starting point to opts.u_max tends to land
%   IPOPT in whichever branch happens to be closest to that large step,
%   which is not necessarily the one that remains feasible all the way
%   down to opts.u_max. A gradual ladder tracks a single branch
%   continuously instead, which is what actually reached the target in
%   practice. If a rung fails, this locally bisects between the last
%   success and that rung before continuing the ladder. If even the
%   "easy" starting point fails to solve, continuation has nothing to
%   warm-start from -- the ORIGINAL error (ME_direct) is surfaced rather
%   than a confusing secondary one.

    u_target = opts.u_max;
    u_easy   = max(2*u_target, u_target + 20);

    opts_step = opts;
    opts_step.u_max = u_easy;
    try
        [~, ~, ~, ~, ~, X_hi, U_hi, Tf_hi] = solve_once(x0, xf, opts_step, params, mode, [], false);
    catch
        rethrow(ME_direct);
    end
    u_hi = u_easy;

    % Linear (not front-loaded) spacing: verified by hand to track the
    % same solution branch all the way to u_target for the T3
    % loop-the-loop case. A denser-near-target ladder was tried and
    % reliably landed IPOPT in a different (faster, ~2.4s, but NOT
    % feasible at the real torque budget) branch instead -- evidence that
    % how the ladder is spaced, not just its endpoints, decides which
    % branch continuation tracks.
    n_rungs = 12;
    ladder = linspace(u_easy, u_target, n_rungs+1);
    ladder = unique([ladder, u_target], 'stable');
    ladder = sort(ladder, 'descend');
    ladder(ladder >= u_hi) = [];   % below the already-solved easy point

    max_local_bisections = 8;
    for i = 1:numel(ladder)
        u_try = ladder(i);
        n_local = 0;
        while true
            opts_step.u_max = u_try;
            warm = struct('X', X_hi, 'U', U_hi, 'Tf', Tf_hi);
            try
                [~, ~, ~, ~, ~, X_try, U_try, Tf_try] = solve_once(x0, xf, opts_step, params, mode, warm, false);
                X_hi = X_try; U_hi = U_try; Tf_hi = Tf_try; u_hi = u_try;
                break;
            catch
                n_local = n_local + 1;
                if n_local > max_local_bisections
                    error('plan_trajectory_casadi:continuation_failed', ...
                        ['plan_trajectory_casadi: direct solve failed and continuation from u_max=%.3f ' ...
                         'stalled between u_max=%.3f (last feasible) and u_max=%.3f rad/s^2 (target %.3f) ' ...
                         'after %d local bisections.'], ...
                        u_easy, u_hi, u_try, u_target, max_local_bisections);
                end
                u_try = 0.5 * (u_hi + u_try);   % bisect toward the last known-feasible point
            end
        end
    end

    opts_step.u_max = u_target;
    warm = struct('X', X_hi, 'U', U_hi, 'Tf', Tf_hi);
    [t_traj, x_traj, u_traj, Tf_opt, sol] = solve_once(x0, xf, opts_step, params, mode, warm);
    fprintf('Continuation succeeded: reached u_max=%.3f (started from %.3f, %d-rung ladder).\n', ...
        u_target, u_easy, numel(ladder));
end


%% ========================================================================
%  HELPERS
%% ========================================================================

function [X_init, U_init, Tf_init] = build_initial_guess(f, x0, xf, K, Tf_min, Tf_max, R_eff, g, u_max, k_roll, n_swings, psi_dot_f)
%BUILD_INITIAL_GUESS  Physically-motivated initial guess for the IPOPT
%                     solve: a resonance-pumping control sinusoid,
%                     forward-simulated through the CasADi dynamics with
%                     RK4, then overwritten on (psi, psi_dot) by a shape
%                     that guides IPOPT toward the correct topological
%                     branch. n_swings selects which shape.
%
%   Inputs
%   ------
%   f      : CasADi Function, xdot = f(x, u)
%   x0, xf : boundary reduced states [theta; theta_dot; psi; psi_dot]
%            [rad; rad/s; rad; rad/s]
%   K      : number of collocation points
%   Tf_min, Tf_max : trajectory duration bounds [s]
%   R_eff  : ball-centre orbit radius for this mode [m]
%   g      : gravitational acceleration [m/s^2]
%   u_max  : control saturation limit [rad/s^2]
%   k_roll : two-rail rolling factor 1 + I_ball/(m*r_roll^2)  [-]
%   n_swings : number of pump-up swings to write into the psi guess.
%            0 = the monotone cosine ramp (see below).
%
%   Outputs
%   -------
%   X_init  : initial state guess (nx x K)
%   U_init  : initial control guess (1 x K) [rad/s^2]
%   Tf_init : initial duration guess [s]
%
%   WHICH SHAPE, AND WHY IT MATTERS. The monotone cosine ramp (n_swings
%   = 0) walks psi from x0(3) to xf(3) without ever reversing. For a
%   manoeuvre the actuator can reach in one go that is exactly the right
%   hint. For a loop-the-loop it is the wrong topological class: the
%   plan that actually exists at the real torque budget swings the ball
%   back and forth several times to build energy before committing, and
%   a guess that never reverses sign starts IPOPT on the far side of
%   that structure. n_swings > 0 writes the pumping in explicitly --
%   n_swings growing half-cycles at the natural frequency, then the
%   revolution.
%
%   The natural frequency is the ball's own pendulum frequency on a HELD
%   hoop. With the ball carried on two O-rings its effective inertia
%   about the hoop axis is m*R_eff^2*k_roll (rolling_matrices.m line 50),
%   against the gravity torque m*g*R_eff*sin(psi), so
%       omega_n = sqrt(g / (k_roll*R_eff)).
%   This function previously hard-coded the 5/7 of a sphere rolling on a
%   single point of contact, i.e. k_roll = 1.4 -- 5.3 % high on this
%   prototype (8.747 against 8.307 rad/s), and wrong for the same reason
%   ball_radius would be wrong in place of r_roll anywhere else.

    nx = numel(x0);

    omega_n = sqrt(g / (k_roll * R_eff));
    T_n     = 2*pi / omega_n;

    if n_swings > 0
        % Enough time for the swings plus the revolution itself, kept
        % inside the caller's bracket. Half a natural period is the
        % order of magnitude of the revolution once the energy is there.
        Tf_init = min(max((n_swings + 0.5) * T_n, Tf_min), Tf_max);
    else
        Tf_init = min(1.5, 0.5*(Tf_min + Tf_max));   % 1.5 s cap is arbitrary, avoids an initial Tf guess large enough to make IPOPT fail
    end

    t_init  = linspace(0, Tf_init, K);
    dt_init = Tf_init / (K-1);

    % Control: resonance pumping sinusoid
    U_init = u_max * sin(2*pi * t_init / T_n);

    % State init: simulate forward with RK4 (for theta, theta_dot, psi_dot).
    % Forward Euler was tried first but gave a much less accurate initial guess.
    X_init = zeros(nx, K);
    X_init(:,1) = x0;
    for k = 1:K-1
        u_k = U_init(k);
        k1 = full(f(X_init(:,k),                 u_k));
        k2 = full(f(X_init(:,k) + 0.5*dt_init*k1, u_k));
        k3 = full(f(X_init(:,k) + 0.5*dt_init*k2, u_k));
        k4 = full(f(X_init(:,k) +     dt_init*k3, u_k));
        X_init(:,k+1) = X_init(:,k) + (dt_init/6) * (k1 + 2*k2 + 2*k3 + k4);
    end

    if nargin >= 12 && ~isempty(psi_dot_f)
        % Cubic Hermite on psi: value and SLOPE prescribed at both ends,
        % zero at t0 (the ball starts at rest) and psi_dot_f at Tf (the
        % release speed N = 0 demands). The cosine ramp below cannot do
        % this -- it is flat at both ends by construction.
        s  = linspace(0, 1, K);
        p0 = x0(3);  p1 = xf(3);
        m0 = 0;      m1 = psi_dot_f * Tf_init;   % d(psi)/ds = Tf * d(psi)/dt
        h00 =  2*s.^3 - 3*s.^2 + 1;   h10 = s.^3 - 2*s.^2 + s;
        h01 = -2*s.^3 + 3*s.^2;       h11 = s.^3 - s.^2;
        psi_guess     = h00*p0 + h10*m0 + h01*p1 + h11*m1;
        d00 =  6*s.^2 - 6*s;          d10 = 3*s.^2 - 4*s + 1;
        d01 = -6*s.^2 + 6*s;          d11 = 3*s.^2 - 2*s;
        psi_dot_guess = (d00*p0 + d10*m0 + d01*p1 + d11*m1) / Tf_init;
    elseif n_swings > 0
        [psi_guess, psi_dot_guess] = pumping_guess(t_init, Tf_init, x0(3), xf(3), omega_n, n_swings);
    else
        % Smooth (cosine-shaped) ramp for psi guides IPOPT toward the correct
        % topological branch (full loop in the -psi direction), with zero slope
        % at both endpoints so the initial guess doesn't kink at t0/Tf.
        tau = linspace(0, 1, K);
        psi_guess     = x0(3) + (xf(3) - x0(3)) * 0.5 * (1 - cos(pi*tau));
        psi_dot_guess = (xf(3) - x0(3)) * 0.5 * (pi/Tf_init) * sin(pi*tau);
    end
    X_init(3, :) = psi_guess;
    X_init(4, :) = psi_dot_guess;
end


function [psi_g, psi_dot_g] = pumping_guess(t, Tf, psi_0, psi_f, omega_n, n_swings)
%PUMPING_GUESS  psi(t) that oscillates at omega_n with linearly growing
%               amplitude for n_swings swings, then runs to psi_f.
%
%   Inputs
%   ------
%   t        : time grid (1 x K)                                     [s]
%   Tf       : trajectory duration                                   [s]
%   psi_0    : initial ball angle                                  [rad]
%   psi_f    : target ball angle (sign sets the loop direction)    [rad]
%   omega_n  : ball pendulum frequency on a held hoop            [rad/s]
%   n_swings : number of pump-up swings                              [-]
%
%   Outputs
%   -------
%   psi_g, psi_dot_g : guess and its exact derivative (1 x K)
%
%   Amplitude grows LINEARLY to 0.8*pi over the pump phase. Linear
%   because that is what constant-energy-per-swing pumping looks like
%   near the small-angle regime, not because it is the optimum; 0.8*pi
%   because the last swing must come close to the top without the guess
%   itself asserting that it gets over -- the revolution is the next
%   segment's job. Both are guess-shaping numbers: they change where
%   IPOPT starts, never what it is allowed to converge to.

    t_pump   = n_swings * (2*pi / omega_n);
    t_pump   = min(t_pump, 0.8 * Tf);        % always leave room for the revolution
    dir      = sign(psi_f - psi_0);
    if dir == 0, dir = -1; end
    amp_max  = 0.8 * pi;

    psi_g     = zeros(size(t));
    psi_dot_g = zeros(size(t));

    is_pump = t <= t_pump;
    tp      = t(is_pump);

    % Growing-amplitude oscillation about psi_0, launched toward psi_f.
    amp     = amp_max * (tp / t_pump);
    amp_dot = amp_max / t_pump;
    psi_g(is_pump)     = psi_0 + dir * amp .* sin(omega_n * tp);
    psi_dot_g(is_pump) = dir * (amp_dot * sin(omega_n * tp) ...
                                + amp .* omega_n .* cos(omega_n * tp));

    % Revolution: cosine ramp from the end of the pump phase to psi_f,
    % zero slope at Tf so the guess meets the terminal condition cleanly.
    is_loop = ~is_pump;
    if any(is_loop)
        tl        = t(is_loop);
        psi_start = psi_g(find(is_pump, 1, 'last'));
        span      = Tf - t_pump;
        s         = (tl - t_pump) / span;
        psi_g(is_loop)     = psi_start + (psi_f - psi_start) * 0.5 .* (1 - cos(pi*s));
        psi_dot_g(is_loop) = (psi_f - psi_start) * 0.5 * (pi/span) .* sin(pi*s);
    end
end


function report_failed_solve(opti, X, U, Tf, K, u_max)
%REPORT_FAILED_SOLVE  Prints and plots the latest (infeasible/non-converged)
%                     iterate after an IPOPT solve failure, for debugging.
%
%   Inputs
%   ------
%   opti  : CasADi Opti stack, in its last-attempted state
%   X, U, Tf : the Opti decision variables (state, control, duration)
%   K     : number of collocation points
%   u_max : control saturation limit, for the plot bounds [rad/s^2]
%
%   Outputs: none (prints diagnostics, creates a figure)

    fprintf('\nSolve failed.\n');
    debug_X  = opti.debug.value(X);
    debug_U  = opti.debug.value(U);
    debug_Tf = opti.debug.value(Tf);
    fprintf('  Tf debug      = %.3f s\n', debug_Tf);
    fprintf('  psi range     = [%.1f, %.1f] deg\n', ...
            rad2deg(min(debug_X(3,:))), rad2deg(max(debug_X(3,:))));
    fprintf('  psi_dot range = [%.1f, %.1f] rad/s\n', ...
            min(debug_X(4,:)), max(debug_X(4,:)));
    fprintf('  theta range   = [%.1f, %.1f] rad\n', ...
            min(debug_X(1,:)), max(debug_X(1,:)));
    fprintf('  u max abs     = %.2f rad/s^2\n', max(abs(debug_U)));

    figure('Name', 'Failed iterate');
    td = linspace(0, debug_Tf, K);
    subplot(4,1,1); plot(td, debug_X(1,:)); ylabel('\theta [rad]'); grid on;
    subplot(4,1,2); plot(td, rad2deg(debug_X(3,:))); ylabel('\psi [deg]'); grid on;
    subplot(4,1,3); plot(td, debug_X(4,:)); ylabel('d\psi/dt [rad/s]'); grid on;
    subplot(4,1,4); plot(td, debug_U); ylabel('u [rad/s^2]'); xlabel('t [s]'); grid on;
    yline(u_max, '--r'); yline(-u_max, '--r');
end
