function plot_energy(params, t_all, X_all, te_all)
%PLOT_ENERGY  Plots kinetic/potential/total energy and the relative
%             energy-conservation error over a run.
%
%   Inputs
%   ------
%   params : physical parameters from ball_hoop_params()
%   t_all  : time vector (N x 1)                                   [s]
%   X_all  : state history (N x 6), [r, r_dot, theta, theta_dot, psi, psi_dot]
%            [m, m/s, rad, rad/s, rad, rad/s]
%   te_all : times of mode-transition events (M x 1)                [s]
%
%   Outputs: none (creates a figure)

    figure('Name', 'Energy Conservation');

    % Compute energy
    [E_total, E_kin, E_pot, E_error] = compute_energy(t_all, X_all, params, []);

    % Subplot 1: energy components
    subplot(2, 1, 1);
    plot(t_all, E_kin, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Kinetic');
    hold on;
    plot(t_all, E_pot, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Potential');
    plot(t_all, E_total, 'k-', 'LineWidth', 2, 'DisplayName', 'Total');

    % Mark mode-transition events
    if ~isempty(te_all)
        for i = 1:length(te_all)
            xline(te_all(i), '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
        end
    end

    grid on;
    xlabel('Time [s]');
    ylabel('Energy [J]');
    title('Energy Components');
    legend('Location', 'best');

    % Subplot 2: relative error
    subplot(2, 1, 2);
    plot(t_all, E_error, 'k-', 'LineWidth', 1.5);
    hold on;

    % Mark mode-transition events
    if ~isempty(te_all)
        for i = 1:length(te_all)
            xline(te_all(i), '--', 'Color', [0.5 0.5 0.5]);
        end
    end

    % Acceptable band (1%)
    yline(1, 'r--', 'LineWidth', 1);
    yline(-1, 'r--', 'LineWidth', 1);

    grid on;
    xlabel('Time [s]');
    ylabel('Energy Error [%]');
    title(sprintf('Energy Conservation Error (max = %.2f%%)', max(abs(E_error))));

    % Console summary
    fprintf('\n--- Energy Validation ---\n');
    fprintf('Initial energy: %.6f J\n', E_total(1));
    fprintf('Final energy:   %.6f J\n', E_total(end));
    fprintf('Max error:      %.4f %%\n', max(abs(E_error)));
    fprintf('Mean error:     %.4f %%\n', mean(abs(E_error)));

    if max(abs(E_error)) < 1
        fprintf('Energy conservation OK (< 1%%)\n');
    elseif max(abs(E_error)) < 5
        fprintf('Energy conservation acceptable (< 5%%)\n');
    else
        fprintf('Energy conservation FAILED (> 5%%)\n');
    end
end