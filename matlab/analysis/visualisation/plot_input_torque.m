function plot_input_torque(t_all, tau_vals)
%PLOT_INPUT_TORQUE  Plots the applied motor torque over time.
%
%   Inputs
%   ------
%   t_all    : time vector (N x 1)                    [s]
%   tau_vals : applied motor torque history (N x 1)    [N.m]
%
%   Outputs: none (creates a figure)

    figure('Name', 'Input Torque');
    plot(t_all, tau_vals, 'LineWidth', 1.5); hold on;
    grid on;
    xlabel('Time [s]');
    ylabel('Torque τ [N·m]');
    title('Applied Motor Torque');

end 