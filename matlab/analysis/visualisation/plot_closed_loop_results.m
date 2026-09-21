function plot_closed_loop_results(run, params)
%PLOT_CLOSED_LOOP_RESULTS  Standard 4-panel visualization of a closed-loop run.
%
% Displays in a single figure:
%   (1) Tracking: psi(t) vs reference psi_ref(t)
%   (2) Tracking error: e(t) = psi_ref(t) - psi(t)
%   (3) Control torque: tau(t) with saturation bounds
%   (4) Phase portrait: psi_dot vs psi, with the reference trajectory
%
% Inputs:
%   run    : struct built by build_run, obeying results/RUN_SCHEMA.md
%   params : physical parameters (used for unit conversions, optional)

    if nargin < 2, params = run.params; end

    % --- Extract for readability ---
    t        = run.t_s;
    r        = run.x(:, 1);
    r_dot    = run.x(:, 2);
    theta    = run.x(:, 3);
    theta_dot= run.x(:, 4);
    psi      = run.x(:, 5);
    psi_dot  = run.x(:, 6);
    psi_ref  = run.ref_psi_rad;
    tau      = run.tau_Nm;
    tau_max  = run.cfg.tau_max;
    % Angles continus pour l'affichage
    psi_plot     = unwrap(psi);
    psi_ref_plot = unwrap(psi_ref);
    
    psi_deg = rad2deg(psi_plot);
    ref_deg = rad2deg(psi_ref_plot);

    % Erreur continue, cohérente avec les angles déroulés
    err_deg = rad2deg(psi_ref_plot - psi_plot);
    
    %% --- Build figure ---
    figure('Name', run.name, 'Position', [0, 0, 800, 800]);
    
    % (1) Tracking
    subplot(2, 2, 1);
    plot(t, ref_deg, 'k--', 'LineWidth', 1.2, 'DisplayName', '\psi_{ref}');
    hold on;
    plot(t, psi_deg, 'b-',  'LineWidth', 1.5, 'DisplayName', '\psi');
    grid on; xlabel('Time [s]'); ylabel('\psi [deg]');
    title('Tracking');
    legend('Location', 'best');
    
    % (2) Error
    subplot(2, 2, 2);
    plot(t, err_deg, 'r-', 'LineWidth', 1.5);
    grid on; xlabel('Time [s]'); ylabel('e = \psi_{ref} - \psi  [deg]');
    title('Tracking error');
    yline(0, 'k:');
    
    % (3) Torque with saturation bounds
    subplot(2, 2, 3);
    plot(t, tau, 'm-', 'LineWidth', 1.5);
    hold on;
    yline( tau_max, 'k--', 'DisplayName', '+\tau_{max}');
    yline(-tau_max, 'k--', 'DisplayName', '-\tau_{max}');
    grid on; xlabel('Time [s]'); ylabel('\tau [N.m]');
    title('Control torque');
    ylim([-1.2*tau_max, 1.2*tau_max]);
    
    % (4) Phase portrait
    subplot(2, 2, 4);
    plot(psi_deg, psi_dot, 'b-', 'LineWidth', 1.2);
    hold on;
    plot(psi_deg(1),   psi_dot(1),   'go', 'MarkerFaceColor', 'g', ...
         'DisplayName', 'start');
    plot(psi_deg(end), psi_dot(end), 'rs', 'MarkerFaceColor', 'r', ...
         'DisplayName', 'end');
    plot(ref_deg(end), 0, 'k*', 'MarkerSize', 10, 'DisplayName', 'target');
    grid on; xlabel('\psi [deg]'); ylabel('d\psi/dt [rad/s]');
    title('Phase portrait');
    legend('Location', 'best');
    
    sgtitle(run.name, 'Interpreter', 'none');

    %plot_states(t, r, r_dot, theta, theta_dot, psi, psi_dot);

end