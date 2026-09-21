function plot_error(t, e, m, ctrl_type)
%PLOT_ERROR  Tracking error e(t) = psi_ref - psi with integral indices
%
% Integral performance indices (Graham & Lathrop, 1953;
% Astrom & Hagglund, 1995, sec. 4.7):
%   IAE  = int |e| dt
%   ISE  = int e^2 dt
%   ITAE = int t*|e| dt   (penalizes late errors more heavily)
%   ITSE = int t*e^2 dt
%
%   Inputs
%   ------
%   t         : time vector (N x 1)                    [s]
%   e         : tracking error e(t) = psi_ref - psi (N x 1) [rad]
%   m         : metrics struct with fields .IAE [rad.s], .ISE [rad^2.s],
%               .ITAE [rad.s^2], .ITSE [rad^2.s^2], .e_ss [rad]
%   ctrl_type : controller name string, e.g. 'PID' | 'LQR' | 'TVLQR'
%
%   Outputs: none (creates a figure)

    figure('Name', sprintf('Tracking error - %s', ctrl_type));
    
    e_deg = rad2deg(e);
    
    % --- Main error trace ---
    plot(t, e_deg, 'r-', 'LineWidth', 1.5, 'DisplayName', 'e(t)');
    hold on;
    yline(0, 'k:', 'HandleVisibility', 'off');
    
    % --- Steady-state error band visualization ---
    % Highlight the window used for e_ss computation (last 10%)
    n = length(t);
    ss_start_idx = max(1, round(0.9*n));
    x_ss = [t(ss_start_idx), t(end), t(end), t(ss_start_idx)];
    y_lim = ylim;
    y_ss = [y_lim(1), y_lim(1), y_lim(2), y_lim(2)];
    fill(x_ss, y_ss, [0.9 0.9 1.0], ...
         'EdgeColor', 'none', 'FaceAlpha', 0.4, ...
         'HandleVisibility', 'off');
    % Replot e on top so the fill doesn't cover it
    plot(t, e_deg, 'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    % --- Title with all integral indices ---
    title_str = sprintf(['%s: IAE=%.4f, ISE=%.5f, ', ...
                         'ITAE=%.4f, ITSE=%.5f'], ...
                        ctrl_type, m.IAE, m.ISE, m.ITAE, m.ITSE);
    title(title_str);
    
    % --- Annotation for steady-state error ---
    text(0.98, 0.05, sprintf('e_{ss} = %.3f deg', rad2deg(m.e_ss)), ...
         'Units', 'normalized', ...
         'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'bottom', ...
         'BackgroundColor', 'w', ...
         'EdgeColor', 'k');
    
    grid on;
    xlabel('Time [s]');
    ylabel('e = \psi_{ref} - \psi  [deg]');
    legend('Location', 'best');
end