function plot_command(t, tau, tau_max, m, ctrl_type)
%PLOT_COMMAND  Control effort visualization with saturation diagnostics
%
% Shows tau(t) against actuator limits, with control energy
% (int tau^2 dt) and saturation ratio as quality indicators.
% A heavily saturated controller violates the linear assumptions
% used in LQR/PID design (Astrom & Murray, 2008, sec. 11.3).
%
%   Inputs
%   ------
%   t         : time vector (N x 1)                         [s]
%   tau       : applied motor torque history (N x 1)        [N.m]
%   tau_max   : torque saturation limit                     [N.m]
%   m         : metrics struct with fields .E_ctrl [N^2.m^2.s],
%               .tau_max_used [N.m], .saturation_ratio (dimensionless, 0-1)
%   ctrl_type : controller name string, e.g. 'PID' | 'LQR' | 'TVLQR'
%
%   Outputs: none (creates a figure)

    figure('Name', sprintf('Control effort - %s', ctrl_type));
    
    % --- Highlight saturated regions in the background ---
    sat_mask = abs(tau) >= 0.99 * tau_max;
    if any(sat_mask)
        y_lim = [-1.2*tau_max, 1.2*tau_max];
        % Find contiguous saturation intervals and fill them
        sat_diff = diff([0; sat_mask(:); 0]);
        starts = find(sat_diff == 1);
        ends   = find(sat_diff == -1) - 1;
        for k = 1:length(starts)
            x_sat = [t(starts(k)), t(ends(k)), ...
                     t(ends(k)), t(starts(k))];
            y_sat = [y_lim(1), y_lim(1), y_lim(2), y_lim(2)];
            fill(x_sat, y_sat, [1.0 0.85 0.85], ...
                 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
                 'HandleVisibility', 'off');
            hold on;
        end
    end
    
    % --- Torque trace ---
    plot(t, tau, 'b-', 'LineWidth', 1.5, 'DisplayName', '\tau');
    hold on;
    
    % --- Saturation limits ---
    yline( tau_max, 'r--', 'LineWidth', 1.2, ...
           'Label', '+\tau_{max}', 'LabelHorizontalAlignment', 'left', ...
           'HandleVisibility', 'off');
    yline(-tau_max, 'r--', 'LineWidth', 1.2, ...
           'Label', '-\tau_{max}', 'LabelHorizontalAlignment', 'left', ...
           'HandleVisibility', 'off');
    
    % --- Zero line ---
    yline(0, 'k:', 'HandleVisibility', 'off');
    
    % --- Title with effort metrics ---
    title_str = sprintf(['%s: \\int\\tau^2 dt = %.4f N^2m^2s, ', ...
                         '|\\tau|_{max} = %.3f N\\cdotm, ', ...
                         'sat = %.1f%%'], ...
                        ctrl_type, m.E_ctrl, m.tau_max_used, ...
                        100*m.saturation_ratio);
    title(title_str);
    
    grid on;
    xlabel('Time [s]');
    ylabel('\tau [N\cdotm]');
    ylim([-1.2*tau_max, 1.2*tau_max]);
    legend('Location', 'best');
end