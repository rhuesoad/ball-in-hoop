function result = check_mode_sequence(expected_sequence, mode_intervals)
%CHECK_MODE_SEQUENCE  Compares a hybrid run's realised mode sequence
%                     against the sequence a scenario declared it expects.
%
%   result = check_mode_sequence(expected_sequence, mode_intervals)
%
%   The single most useful diagnostic for a hybrid system: rather than
%   requiring a human to eyeball "Mode transition: ..." console lines,
%   this reports exactly where the realised trajectory first departs
%   from the declared plan.
%
%   Inputs
%   ------
%   expected_sequence : cell array of mode strings, in order, e.g.
%                        {'rolling_out', 'free_fall', 'rolling_in_inside'}
%   mode_intervals     : struct array from ball_hoop_ode.m's sol.mode_intervals,
%                        fields .t_start, .t_end [s], .mode
%
%   Outputs
%   -------
%   result : struct with field
%              .matches : logical, true iff the realised sequence equals
%                         expected_sequence exactly (same modes, same order,
%                         same count)
%            and, only when .matches is false:
%              .divergence_index : index (1-based) of the first interval
%                                  that disagrees, into whichever of the
%                                  two sequences is shorter at that point
%              .expected_mode    : expected_sequence{divergence_index}, or
%                                  '' if the realised sequence ran out first
%              .actual_mode      : mode_intervals(divergence_index).mode, or
%                                  '' if the expected sequence ran out first
%              .divergence_time_s : mode_intervals(divergence_index).t_start [s],
%                                   or the last realised interval's t_end if
%                                   the realised sequence ran out first
%              .divergence_state  : not populated here (this function only
%                                   sees mode_intervals, not the state
%                                   trajectory) -- run_scenario.m fills this
%                                   in from sol.X at .divergence_time_s

    actual_sequence = {mode_intervals.mode};

    n_expected = numel(expected_sequence);
    n_actual   = numel(actual_sequence);
    n_common   = min(n_expected, n_actual);

    divergence_index = 0;
    for k = 1:n_common
        if ~strcmp(expected_sequence{k}, actual_sequence{k})
            divergence_index = k;
            break;
        end
    end
    if divergence_index == 0 && n_expected ~= n_actual
        divergence_index = n_common + 1;
    end

    if divergence_index == 0
        result.matches = true;
        return;
    end

    result.matches           = false;
    result.divergence_index  = divergence_index;

    if divergence_index <= n_expected
        result.expected_mode = expected_sequence{divergence_index};
    else
        result.expected_mode = '';
    end

    if divergence_index <= n_actual
        result.actual_mode        = actual_sequence{divergence_index};
        result.divergence_time_s  = mode_intervals(divergence_index).t_start;
    else
        result.actual_mode        = '';
        result.divergence_time_s  = mode_intervals(end).t_end;
    end
end
