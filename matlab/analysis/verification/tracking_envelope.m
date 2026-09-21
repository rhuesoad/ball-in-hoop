function env = tracking_envelope(mode, params, cfg, opts)
%TRACKING_ENVELOPE  Largest ball oscillation psi the bench can produce,
%                   as a function of the drive frequency.
%
%   env = tracking_envelope()
%   env = tracking_envelope(mode, params, cfg, opts)
%
%   The question. We drive the hoop with a sinusoid theta(t) and the ball
%   answers with a sinusoid psi(t) at the SAME frequency -- the system is
%   linear for small oscillations, so no other frequency appears. How
%   large can that psi be, at each frequency?
%
%   The answer separates into two independent stages, and keeping them
%   apart is the whole point of this function:
%
%   STAGE 1 -- what the hoop can do. Nothing to do with the ball. A
%   sinusoid theta(t) = A_theta*sin(w*t) simultaneously commits an
%   amplitude A_theta, a speed A_theta*w and an acceleration A_theta*w^2.
%   The bench bounds all three separately, so
%
%       A_theta_max(w) = min( theta_max , theta_dot_max/w , u_max/w^2 )
%
%   The three bind in different places: travel at low frequency, torque at
%   high frequency, speed in between. That is why no single one of them
%   answers the question.
%
%   STAGE 2 -- what the ball does with it. A plain Bode magnitude of
%
%       Psi(s)/Theta(s) = -(M21*s^2 + C21*s) / (M22*s^2 + C22*s + k_g)
%
%   with k_g = m*g*R_eff. Damping is kept, unlike an earlier version of
%   this file that dropped it: without C22 the magnitude diverges at the
%   ball's natural frequency and the most interesting part of the plot is
%   a meaningless infinity. With it the peak is finite, ~1/(2*zeta), and
%   the resonance region becomes usable -- which matters, because that is
%   exactly where the achievable amplitude is largest.
%
%   THE RESULT is the product:  A_psi_max(w) = |H(jw)| * A_theta_max(w)
%
%   Inputs
%   ------
%   mode   : rolling mode string          (default 'rolling_out')
%   params : physical parameters          (default ball_hoop_params())
%   cfg    : numerical settings           (default sim_config())
%   opts   : (optional) struct, any subset of
%       .u_max          hoop acceleration bound [rad/s^2]. Default
%                       cfg.tau_max / M(1,1), i.e. the torque limit
%                       divided by the inertia the motor actually sees
%                       WITH the ball on the track (not
%                       params.hoop_motor_inertia, which is the hoops and
%                       rotor alone).
%                       CAUTION: cfg.tau_max is itself marked PROVISIONAL
%                       in sim_config.m -- it assumes a 17.5 A driver
%                       limit that is not configured yet (8 A today), and
%                       a regression of the E3 currents suggests the real
%                       inertia is ~1.65x M(1,1). Both push u_max DOWN,
%                       by a combined factor of about 5. Pass this
%                       explicitly for anything quotable.
%       .theta_dot_max  hoop speed bound [rad/s]. Default 2 turns/s =
%                       12.566, the ODrive vel_limit the E3 campaign ran
%                       with. A CONFIGURED limit, far below the motor's
%                       electrical ceiling (~679 rad/s no-load at 24 V).
%       .theta_max      hoop travel bound [rad]. Default Inf.
%                       TODO(measure): the real bound is whatever the
%                       wiring to the hoop tolerates before winding up.
%                       Left out rather than guessed -- but if finite it
%                       dominates at low frequency.
%       .f_hz           frequency grid [Hz]
%       .plot           draw the three-panel figure  (default true)
%
%   Outputs
%   -------
%   env : struct with .f_hz, .H_mag (stage 2), .A_theta_max_rad (stage 1),
%         .A_psi_max_deg (the product), the individual stage-1 limits, and
%         .f0_hz, .zeta.

    if nargin < 1 || isempty(mode),   mode   = 'rolling_out';      end
    if nargin < 2 || isempty(params), params = ball_hoop_params(); end
    if nargin < 3 || isempty(cfg),    cfg    = sim_config();       end
    if nargin < 4 || isempty(opts),   opts   = struct();           end

    [M, C, ~] = rolling_matrices(params, 0, mode);
    [~, R_eff] = hoop_geometry(mode, params);
    M21 = M(2,1);  M22 = M(2,2);
    C21 = C(2,1);  C22 = C(2,2);
    k_g = params.ball_mass * params.gravity * R_eff;

    defaults = struct( ...
        'u_max',         cfg.tau_max / M(1,1), ...
        'theta_dot_max', 2 * 2*pi, ...
        'theta_max',     Inf, ...
        'f_hz',          logspace(log10(0.02), log10(5), 600), ...
        'plot',          true);
    fn = fieldnames(defaults);
    for i = 1:numel(fn)
        if ~isfield(opts, fn{i}), opts.(fn{i}) = defaults.(fn{i}); end
    end

    w = 2*pi*opts.f_hz(:);

    % --- Stage 1: what the hoop can do ---
    A_travel = opts.theta_max     * ones(size(w));
    A_speed  = opts.theta_dot_max ./ w;
    A_accel  = opts.u_max         ./ w.^2;
    A_theta  = min([A_travel, A_speed, A_accel], [], 2);

    % --- Stage 2: what the ball does with it ---
    s = 1i * w;
    H = -(M21*s.^2 + C21*s) ./ (M22*s.^2 + C22*s + k_g);
    H_mag = abs(H);

    % --- Result ---
    A_psi = H_mag .* A_theta;

    w0   = sqrt(k_g / M22);
    zeta = C22 / (2 * w0 * M22);

    env = struct( ...
        'mode',            mode, ...
        'f_hz',            opts.f_hz(:), ...
        'H',               H, ...
        'H_mag',           H_mag, ...
        'M',               M, ...
        'C',               C, ...
        'A_theta_max_rad', A_theta, ...
        'A_travel_rad',    A_travel, ...
        'A_speed_rad',     A_speed, ...
        'A_accel_rad',     A_accel, ...
        'A_psi_max_deg',   rad2deg(A_psi), ...
        'f0_hz',           w0/(2*pi), ...
        'zeta',            zeta, ...
        'u_max',           opts.u_max, ...
        'theta_dot_max',   opts.theta_dot_max, ...
        'theta_max',       opts.theta_max);

    if opts.plot
        plot_envelope(env);
    end
end


function plot_envelope(env)
%PLOT_ENVELOPE  The two stages and their product, stacked so the shape of
%               the answer can be traced back to which stage caused it.
    figure('Name', sprintf('Tracking envelope -- %s', env.mode), ...
        'Color', 'w', 'Position', [100 100 760 860]);
    f = env.f_hz;

    % --- Stage 1 ---
    ax1 = subplot(3,1,1); hold(ax1,'on');
    h = []; lab = {};
    if isfinite(env.theta_max)
        h(end+1) = plot(ax1, f, rad2deg(env.A_travel_rad), '--', 'Color',[0.30 0.65 0.30]);
        lab{end+1} = 'travel bound';
    end
    h(end+1) = plot(ax1, f, rad2deg(env.A_speed_rad), '--', 'Color',[0.20 0.45 0.80]);
    lab{end+1} = 'speed bound';
    h(end+1) = plot(ax1, f, rad2deg(env.A_accel_rad), '--', 'Color',[0.85 0.35 0.10]);
    lab{end+1} = 'torque bound';
    h(end+1) = plot(ax1, f, rad2deg(env.A_theta_max_rad), 'k-', 'LineWidth', 2);
    lab{end+1} = 'A_\theta max';
    finish(ax1, env, '\bf Stage 1: what the hoop can do', 'hoop amplitude A_\theta  [deg]', [1 1e5]);
    legend(ax1, h, lab, 'Location','southwest');

    % --- Stage 2 ---
    ax2 = subplot(3,1,2);
    plot(ax2, f, env.H_mag, 'k-', 'LineWidth', 2);
    finish(ax2, env, ...
        sprintf('\\bf Stage 2: what the ball does with it   (\\zeta = %.3f, peak \\approx %.0f)', ...
                env.zeta, 1/(2*env.zeta)), ...
        '|\psi/\theta|   [-]', []);

    % --- Product ---
    ax3 = subplot(3,1,3);
    plot(ax3, f, env.A_psi_max_deg, 'k-', 'LineWidth', 2);
    finish(ax3, env, '\bf Result: largest achievable ball oscillation', ...
        'A_\psi max  [deg]', [0.1 1e3]);
    xlabel(ax3, 'drive frequency  [Hz]');
end


function finish(ax, env, ttl, ylab, ylim_)
%FINISH  Shared log-log formatting, with the ball's natural frequency marked.
    set(ax, 'XScale','log', 'YScale','log');
    grid(ax, 'on');
    xlim(ax, [env.f_hz(1), env.f_hz(end)]);
    if ~isempty(ylim_), ylim(ax, ylim_); end
    xline(ax, env.f0_hz, ':', sprintf('f_0 = %.2f Hz', env.f0_hz), ...
        'LineWidth', 1.2, 'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
    ylabel(ax, ylab);
    title(ax, ttl);
end
