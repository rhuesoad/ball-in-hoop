function [new_mode, x0] = transition_state(current_mode, x_pre, params, cfg)
%TRANSITION_STATE  Determines the next dynamic mode after an event and
%                  prepares the initial state for the next segment.
%
%   [new_mode, x0] = transition_state(current_mode, x_pre, params, cfg)
%
%   Inputs
%   ------
%   current_mode : mode string at the time of the event
%   x_pre        : state vector at the event [r; r_dot; theta; theta_dot; psi; psi_dot]
%   params       : physical parameters struct (see ball_hoop_params)
%   cfg          : numerical settings struct (see sim_config)
%
%   Outputs
%   -------
%   new_mode : mode string for the next segment
%   x0       : initial state for the next segment (nudged to avoid re-detection)

% Every path out of a rolling mode goes to free_fall, and the ball's spin
% has to survive that flight: it is locked to the other states while
% rolling but frozen and independent once airborne, so it is computed
% here, on the last state that still satisfies the rolling constraint,
% and handed to handle_impact.m at the far end. See ball_spin_store.m.
if ismember(current_mode, {'rolling_out', 'rolling_in_outside', 'rolling_in_inside'})
    ball_spin_store('set', ball_spin_from_rolling(x_pre, params, current_mode));
end

switch current_mode
    case 'rolling_out'
        new_mode = 'free_fall';
        x0    = x_pre;
        x0(1) = x0(1) - cfg.TOL_NUDGE;  % move slightly inward off the outer hoop
        fprintf('Mode transition: %s -> %s\n', current_mode, new_mode);

    case 'free_fall'
        [new_mode, x0] = handle_impact(x_pre, params, cfg);
        if ~strcmp(new_mode, 'free_fall')
            fprintf('Mode transition: %s -> %s\n', current_mode, new_mode);
        end

    case 'rolling_in_outside'
        new_mode = 'free_fall';
        x0    = x_pre;
        x0(1) = x0(1) + cfg.TOL_NUDGE;  % move slightly outward off the inner hoop
        fprintf('Mode transition: %s -> %s\n', current_mode, new_mode);

    case 'rolling_in_inside'
        new_mode = 'free_fall';
        x0    = x_pre;
        x0(1) = x0(1) - cfg.TOL_NUDGE;  % move slightly inward off the inner hoop
        fprintf('Mode transition: %s -> %s\n', current_mode, new_mode);

    otherwise
        error('transition_state: unknown mode ''%s''.', current_mode);
end
end
