function out = fig1_bode_margins(varargin)
%FIG1_BODE_MARGINS  Defence slide: designed loop vs implemented loop.
%
%   out = fig1_bode_margins()
%   out = fig1_bode_margins('K', [1 4 18 5.5], 'save_dir', d)
%
%   Draws the open-loop Bode of the two loops loop_margins.m builds -- the
%   continuous state feedback the LQR was designed against, and the same
%   feedback seen through the 50 Hz sampler, the 19.3 ms camera transport
%   and the psi_dot difference-and-filter chain the bench actually runs.
%   The gain crossover and the phase margin are marked on the trace, and
%   the gain margin is marked at the phase crossover.
%
%   The magnitude and phase panels are built by hand from squeeze(bode(...))
%   rather than by bodeplot, so the slide keeps one style with the rest of
%   the report: Hz on the abscissa (bodeplot insists on rad/s unless the
%   toolbox preference is changed globally), no auto title, no boxed axes,
%   and the two curves distinguished by colour and dash rather than by a
%   legend the audience has to decode.
%
%   TYPOGRAPHY. The report figures are drawn by matplotlib with usetex, so
%   every label is set by LaTeX in Computer Modern. Here the same is done by
%   switching the interpreter rather than by naming a font: MATLAB silently
%   falls back to Helvetica when given a font it cannot resolve, which is
%   what a bare FontName would do. Note that MATLAB's LaTeX interpreter runs
%   a minimal preamble -- siunitx is NOT available, so units are written out
%   in plain brackets.
%
%   ON THE GAIN SIGN. K is exposed because the two sides of this project
%   count psi in opposite directions (loop_margins.m docstring, and
%   bench/t1_rebroussement_scenarii.m's to_matlab). The default below is
%   the MATLAB-convention gain. Negating the psi and psi_dot entries gives
%   a loop whose DESIGNED form already has a negative gain margin -- i.e.
%   unstable before any delay is added -- which is the signature of having
%   applied that conversion once too often.
%
%   Inputs (name/value, all optional)
%   ------
%   'K'        : 1x4 state feedback [theta, theta_dot, psi, psi_dot].
%                Default [1 4 18 5.5].
%   'Ts'       : control period [s]. Default 1/50.
%   'tau'      : camera transport delay [s]. Default 0.0193.
%   'fc'       : psi_dot low-pass cutoff [Hz]. Default 12.
%   'mode'     : rolling surface. Default 'rolling_out'.
%   'f_lim'    : [f_min f_max] drawn, in Hz. Default [0.05 50].
%   'save_dir' : output directory. Default results/figures/defence.
%
%   Outputs
%   -------
%   out : the loop_margins.m result struct, plus .pdf_path and .png_path.

    p = inputParser;
    p.addParameter('K',    [1, 4, 18, 5.5], @(v) numel(v) == 4);
    p.addParameter('Ts',   1/50, @isscalar);
    p.addParameter('tau',  0.0193, @isscalar);
    p.addParameter('fc',   12, @isscalar);
    p.addParameter('mode', 'rolling_out', @ischar);
    p.addParameter('f_lim', [0.05 50], @(v) numel(v) == 2);
    p.addParameter('save_dir', '', @(s) ischar(s) || isstring(s));
    p.parse(varargin{:});
    opt = p.Results;

    project_root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    if isempty(opt.save_dir)
        opt.save_dir = fullfile(project_root, 'results', 'figures', 'defence');
    end
    if ~exist(opt.save_dir, 'dir'), mkdir(opt.save_dir); end

    %% --- Style ---------------------------------------------------------
    C_BLUE = [0.00 0.29 0.58];    % ULB blue -- implemented loop
    C_GREY = [0.55 0.58 0.62];    % designed loop
    C_RED  = [0.69 0.23 0.18];    % limits and margins
    FSIZE  = 11;

    % LaTeX everywhere, restored on exit -- including on error, so the
    % setting never leaks into the next figure drawn in this session.
    set(groot, 'defaultTextInterpreter',          'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter',        'latex');
    restore_interpreter = onCleanup(@() set(groot, ...
        'defaultTextInterpreter',          'tex', ...
        'defaultAxesTickLabelInterpreter', 'tex', ...
        'defaultLegendInterpreter',        'tex'));

    %% --- The two loops -------------------------------------------------
    % loop_margins.m is the single derivation of both transfer functions
    % and of the margins quoted below; nothing is recomputed here.
    out = loop_margins('K', opt.K, 'Ts', opt.Ts, 'tau', opt.tau, ...
                       'fc', opt.fc, 'mode', opt.mode, 'plot', false);

    w = logspace(log10(2*pi*opt.f_lim(1)), log10(2*pi*opt.f_lim(2)), 3000);
    f = w / (2*pi);

    [mag_i, ph_i] = bode(out.L_ideal, w);
    [mag_c, ph_c] = bode(out.L,       w);
    mag_i = 20*log10(squeeze(mag_i));   ph_i = squeeze(ph_i);
    mag_c = 20*log10(squeeze(mag_c));   ph_c = squeeze(ph_c);

    f_cp = out.Wcp / (2*pi);            % gain crossover, |L| = 0 dB
    f_cg = out.Wcg / (2*pi);            % phase crossover, arg L = -180 deg
    ph_at_cp  = interp1(f, ph_c,  f_cp);
    mag_at_cg = interp1(f, mag_c, f_cg);

    %% --- Figure, 16:9 --------------------------------------------------
    W_CM = 25.4;  H_CM = 14.3;
    fig = figure('Color', 'w', 'Units', 'centimeters', ...
                 'Position', [2 2 W_CM H_CM], 'Name', 'Loop margins');

    tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    % ---- Magnitude -----------------------------------------------------
    ax1 = nexttile(tl); hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'off');
    set(ax1, 'XScale', 'log', 'FontSize', FSIZE);
    yline(ax1, 0, '-', 'Color', [0.80 0.80 0.80], 'LineWidth', 0.8);
    plot(ax1, f, mag_i, '--', 'Color', C_GREY, 'LineWidth', 1.4);
    plot(ax1, f, mag_c, '-',  'Color', C_BLUE, 'LineWidth', 1.8);

    % Gain margin: the drop from 0 dB down to |L| at the phase crossover.
    % The designed loop has no phase crossover at all, so this segment
    % exists only on the blue curve -- which is the point of the slide.
    if isfinite(f_cg) && f_cg > opt.f_lim(1) && f_cg < opt.f_lim(2)
        plot(ax1, [f_cg f_cg], [mag_at_cg 0], '-', ...
             'Color', C_RED, 'LineWidth', 1.6);
        plot(ax1, f_cg, mag_at_cg, 'o', 'Color', C_RED, ...
             'MarkerFaceColor', C_RED, 'MarkerSize', 5);
        text(ax1, f_cg*1.15, mag_at_cg/2, ...
             sprintf('$\\mathrm{GM} = %.1f$ dB', out.Gm_dB), ...
             'Color', C_RED, 'FontSize', FSIZE-1, ...
             'VerticalAlignment', 'middle');
    end

    % Gain crossover, carried down to the phase panel by the same line.
    xline(ax1, f_cp, ':', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.2);
    plot(ax1, f_cp, 0, 'o', 'Color', C_BLUE, 'MarkerFaceColor', 'w', ...
         'LineWidth', 1.4, 'MarkerSize', 6);

    ylabel(ax1, '$|L|$ [dB]', 'FontSize', FSIZE);
    xlim(ax1, opt.f_lim);
    set(ax1, 'XTickLabel', []);

    % Hand-built key: two colour samples and their captions. A legend box
    % would sit on the curves at this aspect ratio.
    yl  = ylim(ax1);
    y0  = yl(1) + 0.10*diff(yl);
    dy  = 0.075*diff(yl);
    x0  = opt.f_lim(1)*1.10;
    x1  = x0*1.35;
    plot(ax1, [x0 x1], [y0+dy y0+dy], '--', 'Color', C_GREY, 'LineWidth', 1.4);
    plot(ax1, [x0 x1], [y0    y0   ], '-',  'Color', C_BLUE, 'LineWidth', 1.8);
    text(ax1, x1*1.15, y0+dy, 'designed (no delay, no filter)', ...
         'Color', C_GREY, 'FontSize', FSIZE-1, 'VerticalAlignment', 'middle');
    text(ax1, x1*1.15, y0, 'implemented (50 Hz, camera, filter)', ...
         'Color', C_BLUE, 'FontSize', FSIZE-1, 'VerticalAlignment', 'middle');

    % ---- Phase ---------------------------------------------------------
    ax2 = nexttile(tl); hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'off');
    set(ax2, 'XScale', 'log', 'FontSize', FSIZE);
    yline(ax2, -180, '-', 'Color', [0.80 0.80 0.80], 'LineWidth', 0.8);
    plot(ax2, f, ph_i, '--', 'Color', C_GREY, 'LineWidth', 1.4);
    plot(ax2, f, ph_c, '-',  'Color', C_BLUE, 'LineWidth', 1.8);

    % Phase margin: the rise from -180 deg up to arg L at the gain crossover.
    xline(ax2, f_cp, ':', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.2);
    plot(ax2, [f_cp f_cp], [-180 ph_at_cp], '-', ...
         'Color', C_RED, 'LineWidth', 1.6);
    plot(ax2, f_cp, ph_at_cp, 'o', 'Color', C_RED, ...
         'MarkerFaceColor', C_RED, 'MarkerSize', 5);
    text(ax2, f_cp*1.15, 0.5*(-180 + ph_at_cp), ...
         sprintf('$\\mathrm{PM} = %.1f^{\\circ}$', out.Pm_deg), ...
         'Color', C_RED, 'FontSize', FSIZE-1, 'VerticalAlignment', 'middle');
    text(ax2, f_cp, ph_at_cp, sprintf('  $%.2f$ Hz', f_cp), ...
         'Color', [0.25 0.25 0.25], 'FontSize', FSIZE-1, ...
         'VerticalAlignment', 'bottom');

    xlabel(ax2, '$f$ [Hz]', 'FontSize', FSIZE);
    ylabel(ax2, '$\angle L$ [deg]', 'FontSize', FSIZE);
    xlim(ax2, opt.f_lim);
    ylim(ax2, [max(-540, min(ph_c)) - 20, 10]);

    %% --- Export --------------------------------------------------------
    stem = fullfile(opt.save_dir, 'fig1_bode_margins');
    set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [W_CM H_CM], ...
             'PaperPosition', [0 0 W_CM H_CM]);
    print(fig, '-dpdf', '-vector', [stem '.pdf']);
    print(fig, '-dpng', '-r300',   [stem '.png']);
    out.pdf_path = [stem '.pdf'];
    out.png_path = [stem '.png'];

    %% --- Report --------------------------------------------------------
    fprintf('\n=== Figure 1: designed vs implemented loop ===\n');
    fprintf('  K                  : [%+.4f %+.4f %+.4f %+.4f]\n', opt.K);
    fprintf('  Ts / tau / fc      : %.1f ms / %.1f ms / %.1f Hz\n', ...
            1e3*opt.Ts, 1e3*opt.tau, opt.fc);
    fprintf('  gain margin  [dB]  : %8.2f (designed) -> %8.2f (implemented)\n', ...
            out.ideal.Gm_dB, out.Gm_dB);
    fprintf('  phase margin [deg] : %8.2f (designed) -> %8.2f (implemented)\n', ...
            out.ideal.Pm_deg, out.Pm_deg);
    fprintf('  gain crossover [Hz]: %8.3f (designed) -> %8.3f (implemented)\n', ...
            out.ideal.Wcp/(2*pi), f_cp);
    fprintf('  phase crossover[Hz]: %8.3f (designed) -> %8.3f (implemented)\n', ...
            out.ideal.Wcg/(2*pi), f_cg);
    fprintf('  delay margin  [ms] : %8.1f on top of the %.1f already modelled\n', ...
            1e3*out.delay_margin_s, 1e3*(opt.tau + opt.Ts));
    fprintf('  PDF : %s\n', out.pdf_path);
    fprintf('  PNG : %s\n', out.png_path);
end