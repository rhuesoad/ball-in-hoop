function plot_tracking(t, X, psi_ref, m, ctrl_type)
%PLOT_TRACKING  Reference tracking visualization with time-domain metrics
%
% Displays psi(t) and psi_ref(t) overlaid, with annotations for the
% standard step-response indicators (Ogata, 2010, sec. 5-3):
%   - rise time (10%-90%)
%   - peak time and overshoot
%   - settling time (2% band)
%   - steady-state error
%
%   Inputs
%   ------
%   t         : time vector (N x 1)                                [s]
%   X         : state history (N x 6), [r, r_dot, theta, theta_dot, psi, psi_dot]
%               [m, m/s, rad, rad/s, rad, rad/s]
%   psi_ref   : reference ball angle history (N x 1)                [rad]
%   m         : metrics struct with fields .t_rise [s], .t_peak [s],
%               .M_p (percent overshoot, dimensionless), .t_settle [s],
%               .e_ss [rad]
%   ctrl_type : controller name string, e.g. 'PID' | 'LQR' | 'TVLQR'
%
%   Outputs: none (creates a figure)

    figure('Name', sprintf('Tracking - %s', ctrl_type));
    
    psi = X(:, 5);
    psi_deg     = rad2deg(psi);
    psi_ref_deg = rad2deg(psi_ref);
    
    % --- Reference and response ---
    plot(t, psi_ref_deg, 'k--', 'LineWidth', 1.5, ...
         'DisplayName', '\psi_{ref}'); hold on;
    plot(t, psi_deg, 'b-', 'LineWidth', 1.5, ...
         'DisplayName', '\psi');
    
    % --- 2% settling band around the final reference value ---
    % We draw the band only if the reference is approximately constant
    % at the end (otherwise the band is meaningless for tracking).
    psi_target_deg = psi_ref_deg(end);
    psi_init_deg   = psi_deg(1);
    span_deg = abs(psi_target_deg - psi_init_deg);
    
    if span_deg > 1e-6  % meaningful step exists
        band = 0.02 * span_deg;
        yline(psi_target_deg + band, ':', 'Color', [0.5 0.5 0.5], ...
              'HandleVisibility', 'off');
        yline(psi_target_deg - band, ':', 'Color', [0.5 0.5 0.5], ...
              'HandleVisibility', 'off');
    end
    
    % --- Markers for t_r, t_p, t_s ---
    % Only plot markers if the metrics were successfully computed
    if ~isnan(m.t_peak)
        % Find psi value at t_peak
        psi_at_peak = interp1(t, psi_deg, m.t_peak);
        plot(m.t_peak, psi_at_peak, 'rv', 'MarkerSize', 8, ...
             'MarkerFaceColor', 'r', 'DisplayName', ...
             sprintf('Peak (M_p=%.1f%%)', m.M_p));
    end
    
    if ~isnan(m.t_settle)
        xline(m.t_settle, ':', 'Color', [0 0.5 0], 'LineWidth', 1.2, ...
              'Label', sprintf('t_s=%.2fs', m.t_settle), ...
              'LabelVerticalAlignment', 'bottom', ...
              'HandleVisibility', 'off');
    end
    
    % --- Title with key metrics ---
    title_str = sprintf('%s: t_r=%.3fs, M_p=%.1f%%, t_s=%.3fs, e_{ss}=%.2f deg', ...
                       ctrl_type, m.t_rise, m.M_p, m.t_settle, ...
                       rad2deg(m.e_ss));
    title(title_str);
    
    grid on;
    xlabel('Time [s]');
    ylabel('\psi [deg]');
    legend('Location', 'best');
end