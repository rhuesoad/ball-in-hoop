function plot_phase_portrait(X, psi_target, ctrl_type)
%PLOT_PHASE_PORTRAIT  Trajectory in the (psi, psi_dot) plane
%
% The phase portrait reveals qualitative properties invisible in
% time-domain plots: spiral inwards (damped), limit cycles,
% multiple equilibria. Standard tool for nonlinear systems
% (Khalil, 2002, ch. 2).
%
%   Inputs
%   ------
%   X          : state history (N x 6), [r, r_dot, theta, theta_dot, psi, psi_dot]
%                [m, m/s, rad, rad/s, rad, rad/s]
%   psi_target : reference ball angle [rad]
%   ctrl_type  : controller name string, e.g. 'PID' | 'LQR' | 'TVLQR'
%
%   Outputs: none (creates a figure)

    figure('Name', sprintf('Phase portrait - %s', ctrl_type));
    
    psi_deg   = rad2deg(X(:, 5));
    psi_d_deg = rad2deg(X(:, 6));
    
    % --- Trajectory with color gradient (time progression) ---
    % Using a colormap shows direction of motion without arrows,
    % which can clutter spiral trajectories.
    n = length(psi_deg);
    cmap = parula(n);
    
    % Plot as line segments to enable color gradient
    for k = 1:(n-1)
        plot(psi_deg(k:k+1), psi_d_deg(k:k+1), '-', ...
             'Color', cmap(k, :), 'LineWidth', 1.2, ...
             'HandleVisibility', 'off');
        hold on;
    end
    
    % Re-plot a thin overlay for the legend entry
    plot(NaN, NaN, 'b-', 'LineWidth', 1.5, ...
         'DisplayName', 'Trajectory (color = time)');
    
    % --- Start and target markers ---
    plot(psi_deg(1), psi_d_deg(1), 'go', ...
         'MarkerSize', 10, 'LineWidth', 2, ...
         'MarkerFaceColor', 'g', ...
         'DisplayName', 'Start');
    plot(rad2deg(psi_target), 0, 'rp', ...
         'MarkerSize', 14, 'LineWidth', 2, ...
         'MarkerFaceColor', 'r', ...
         'DisplayName', 'Target');
    
    % --- End marker (where the trajectory actually settled) ---
    plot(psi_deg(end), psi_d_deg(end), 'ks', ...
         'MarkerSize', 10, 'LineWidth', 2, ...
         'DisplayName', 'End');
    
    % --- Reference axes through the target ---
    xline(rad2deg(psi_target), ':', 'Color', [0.5 0.5 0.5], ...
          'HandleVisibility', 'off');
    yline(0, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    
    % --- Colorbar to indicate time progression ---
    colormap(parula);
    cb = colorbar;
    cb.Label.String = 'Time progression';
    cb.Ticks = [0, 1];
    cb.TickLabels = {'t_0', 't_{end}'};
    
    grid on;
    xlabel('\psi [deg]');
    ylabel('d\psi/dt [deg/s]');
    title(sprintf('Phase portrait - %s', ctrl_type));
    legend('Location', 'best');
    axis tight;
end