function out = ball_spin_store(action, value)
%BALL_SPIN_STORE  Carries the ball's own spin across a free-flight segment.
%
%   ball_spin_store('reset')          clears it (call once per simulation)
%   ball_spin_store('set', phi_dot)   store the spin at liftoff  [rad/s]
%   phi_dot = ball_spin_store('get')  retrieve it at impact      [rad/s]
%
%   WHY THIS EXISTS. The ball has two distinct rotations: psi_dot, how
%   fast it travels around the hoop, and phi_dot, how fast it turns on
%   its own axis. While it ROLLS the two are locked together by the
%   no-slip constraint, so phi_dot is not an independent state and the
%   6-state vector [r; r_dot; theta; theta_dot; psi; psi_dot] loses
%   nothing by omitting it.
%
%   In FLIGHT that stops being true. Nothing torques the ball about its
%   own centre -- gravity acts at the centre, and the drag is 1e-4 of the
%   weight (analysis/studies/T4/t4c_phase1_liftoff.m) -- so phi_dot is FROZEN at
%   its liftoff value while psi_dot and theta_dot go on changing. The
%   rolling relation therefore cannot be used to recover it at the
%   moment of impact, and the state vector has no room for it.
%
%   It is not integrated, only remembered, so a store is enough and no
%   seventh state is needed. That is the whole trick.
%
%   AND IT MATTERS. On the T4c release at 125 deg the liftoff spin is
%   -192.6 rad/s, and in handle_impact.m's law its term I_b*r_roll*phi_dot
%   is about -8.6e-6 against m*r_roll^2*v- of about 7.6e-6 for a 1 m/s
%   arrival: the SAME ORDER. Setting it to zero, as this codebase did
%   before, is not a small approximation.
%
%   PERSISTENCE. A persistent variable is deliberate: the alternative is
%   threading phi_dot through ball_hoop_ode, detect_events,
%   transition_state and every caller of run_scenario, for a quantity
%   that is constant over exactly one segment. The cost is that it
%   survives between simulations, so run_scenario.m resets it before
%   every run; anything else driving the hybrid integrator directly must
%   do the same.

    persistent phi_dot_stored

    switch action
        case 'reset'
            phi_dot_stored = [];
            out = [];

        case 'set'
            phi_dot_stored = value;
            out = [];

        case 'get'
            if isempty(phi_dot_stored)
                % No liftoff was recorded before this impact -- a
                % simulation started mid-flight, for instance. Zero is
                % the old behaviour, and it is announced rather than
                % applied silently.
                warning('ball_spin_store:no_liftoff_recorded', ...
                    ['No ball spin was stored before this impact; assuming 0 rad/s. ' ...
                     'A scenario that starts in free_fall should set it explicitly.']);
                out = 0;
            else
                out = phi_dot_stored;
            end

        otherwise
            error('ball_spin_store:action', ...
                'action must be ''reset'', ''set'' or ''get'', got ''%s''.', action);
    end
end
