function plot_normal_force(params, t_all, r, psi, psi_d, te_all, Xe_all)
%PLOT_NORMAL_FORCE  Normal contact force on the ball, concave-contact
%                    convention (valid for rolling_out / rolling_in_inside;
%                    see docs/AUDIT.md item 1.2 for the convex-contact sign).
%
%   N = m*(g*cos(psi) + R_eff*psi_dot^2), from the radial free-body
%   equation of a ball on the concave (inward-pushing) side of a track.
%
%   Inputs
%   ------
%   params : physical parameters from ball_hoop_params()
%   t_all  : time vector (N x 1)                                [s]
%   r      : ball radial position history (N x 1)                [m]
%   psi    : ball angle history (N x 1)                          [rad]
%   psi_d  : ball angular velocity history (N x 1)                [rad/s]
%   te_all : times of mode-transition events (M x 1)              [s]
%   Xe_all : states at mode-transition events (M x 6),
%            [r, r_dot, theta, theta_dot, psi, psi_dot]
%
%   Outputs: none (creates a figure)

 % --- 4. Normal Force Validation ---
    figure('Name', 'Normal Force');
    N = params.ball_mass * (params.gravity * cos(psi) + r .* psi_d.^2);

    plot(t_all, N, 'LineWidth', 1.5);
    hold on;

    if ~isempty(Xe_all)
        N_events = params.ball_mass * (params.gravity * cos(Xe_all(:, 5)) + ...
                               Xe_all(:, 1).* Xe_all(:, 6).^2);
        plot(te_all, N_events, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    end

    grid on;
    xlabel('Time [s]');
    ylabel('Normal Force [N]');
    title('Normal Force on Ball (Outer Hoop)');
    legend('Normal force', 'Events', 'Location', 'best');

end 