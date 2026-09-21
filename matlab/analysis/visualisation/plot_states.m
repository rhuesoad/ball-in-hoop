function plot_states(t_all, r, r_d, theta, theta_d, psi, psi_d)
%PLOT_STATES  Plots the six physical state components over time.
%
%   Inputs
%   ------
%   t_all   : time vector (N x 1)                 [s]
%   r       : ball radial position (N x 1)         [m]
%   r_d     : ball radial velocity (N x 1)          [m/s]
%   theta   : hoop angle (N x 1)                    [rad]
%   theta_d : hoop angular velocity (N x 1)          [rad/s]
%   psi     : ball angle (N x 1)                     [rad]
%   psi_d   : ball angular velocity (N x 1)           [rad/s]
%
%   Outputs: none (creates a figure)

% --- 2. State Variables ---
    figure('Name', 'System States');
    
    subplot(3, 2, 1);
    plot(t_all, r, 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('r [m]');
    title('Ball Radial Position');

    subplot(3, 2, 2);
    plot(t_all, r_d, 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('ṙ [m/s]');
    title('Ball Radial Velocity');

    subplot(3, 2, 3);
    plot(t_all, theta, 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('θ [rad]');
    title('Hoop Angle');

    subplot(3, 2, 4);
    plot(t_all, theta_d, 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('θ̇ [rad/s]');
    title('Hoop Angular Velocity');

    subplot(3, 2, 5);
    plot(t_all, psi, 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('ψ [rad]');
    title('Ball Angular Position');

    subplot(3, 2, 6);
    plot(t_all, psi_d, 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('ψ̇ [rad/s]');
    title('Ball Angular Velocity');
end