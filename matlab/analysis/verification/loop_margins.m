function out = loop_margins(varargin)
%LOOP_MARGINS  Gain and phase margins of the loop as IMPLEMENTED, delays
%              and psi_dot filter included.
%
%   out = loop_margins()
%   out = loop_margins('mode', 'rolling_in_inside')
%   out = loop_margins('K', [1.0066 6.1964 87.8026 19.8888])
%
%   The state feedback is designed on the continuous plant alone
%   (linearize_system.m -> lqr_design.m), which has infinite gain margin
%   and 60 deg of phase margin by construction. The bench does not run
%   that loop. Between the plant and the gain sit a 50 Hz sampled
%   controller, a camera whose frame reaches the controller 19.3 ms after
%   the instant it shows, and a psi_dot obtained by differencing two
%   camera frames and low-passing the result. Each is pure phase lag, and
%   phase lag is what the LQR guarantee does not cover.
%
%   This function closes that gap: it builds the loop the bench actually
%   implements and reads the margins off it.
%
%   WHAT IS MODELLED, per feedback channel
%   --------------------------------------
%   theta, theta_dot : unity. They come from the ODrive encoder at the
%                      loop rate, not from the camera. See the CAVEAT
%                      below -- on the v3 bench this is optimistic.
%   psi              : exp(-tau*s), camera transport delay.
%   psi_dot          : exp(-tau*s) * exp(-s*Ts/2) / (1 + s/wc).
%                      Camera transport, then the half-sample lag of the
%                      backward difference (psi_k - psi_{k-1})/dt, then
%                      the first-order low-pass at fc
%                      (bench_common_v3.py, RateEstimator).
%
%   The zero-order hold on the command is exp(-s*Ts/2), the usual
%   half-sample approximation, applied once at the plant input.
%
%   CAVEAT, and it is the important one. theta and theta_dot are taken as
%   perfect here because the encoder is fast. On the v3 bench with
%   USE_MEASURED_THETA = False the controller feeds back its OWN
%   integrator rather than the encoder, and the hoop delivers only 37-71%
%   of the commanded velocity, so the true theta channel has a gain error
%   this function does not represent. The margins below are therefore an
%   upper bound on what the bench really has.
%
%   Inputs (name/value)
%   -------------------
%   'mode'    : 'rolling_out' | 'rolling_in_outside' | 'rolling_in_inside'
%               Default 'rolling_out'.
%   'psi_lin' : linearisation point [rad]. Default 0, or pi for
%               'rolling_in_outside' (the inverted equilibrium).
%   'K'       : 1x4 state-feedback gain [theta, theta_dot, psi, psi_dot],
%               MATLAB sign convention. Default: designed here by
%               lqr_design.m with the weights below.
%               A gain copied from the bench config must have its psi and
%               psi_dot entries negated first -- the two sides count psi
%               in opposite directions (see bench/t1_rebroussement_scenarii.m,
%               to_matlab).
%   'Q', 'R'  : LQR weights used when K is not supplied.
%   'Ts'      : control period [s]. Default 1/50 (LOOP_HZ).
%   'tau'     : camera transport delay [s]. Default 0.0193 (CAM_LATENCY_S).
%   'fc'      : psi_dot low-pass cutoff [Hz]. Default 12 (PSIDOT_FILTER_HZ).
%   'w'       : frequency grid [rad/s]. Default logspace(-1, 3, 4000).
%   'plot'    : draw the Bode of both loops. Default true.
%
%   Outputs
%   -------
%   out : struct with
%           .Gm_dB, .Pm_deg, .Wcg, .Wcp   implemented loop
%           .ideal                        same four, delays/filter removed
%           .delay_margin_s               extra latency the loop tolerates
%                                         before Pm reaches zero
%           .L, .L_ideal                  the loop transfer functions
%           .K, .A, .B                    what was used
%
%   Requires the Control System Toolbox.

    if isempty(ver('control'))
        error('loop_margins:toolbox', ...
            'This function needs the Control System Toolbox (tf, ss, margin, frd).');
    end

    p = inputParser;
    p.addParameter('mode', 'rolling_out', @ischar);
    p.addParameter('psi_lin', [], @(v) isempty(v) || isscalar(v));
    p.addParameter('K', [1, 4, 18, 5.5]);
    p.addParameter('Ts', 1/50, @isscalar);
    p.addParameter('tau', 0.0193, @isscalar);
    p.addParameter('fc', 12, @isscalar);
    p.addParameter('w', logspace(-1, 3, 4000), @isnumeric);
    p.addParameter('plot', true, @(v) islogical(v) || isnumeric(v));
    p.parse(varargin{:});
    opt = p.Results;

    if isempty(opt.psi_lin)
        % The inverted equilibrium is the only one not at the bottom.
        if strcmp(opt.mode, 'rolling_in_outside'), opt.psi_lin = pi;
        else,                                     opt.psi_lin = 0;
        end
    end

    params = ball_hoop_params();
    [A, B] = linearize_system(opt.psi_lin, params, opt.mode);

    if isempty(opt.K)
        scn.state = opt.mode;
        scn.ctrl_params = struct('psi_lin', opt.psi_lin, 'Q', opt.Q, 'R', opt.R);
        ctrl = lqr_design(scn, params, sim_config());
        K = ctrl.K(:).';
    else
        K = opt.K(:).';
    end

    %% --- The two loops ---
    s  = tf('s');
    wc = 2*pi*opt.fc;

    % Four outputs, one input: the loop is closed on the full state, so
    % every state must come out of the plant to be delayed individually.
    G = ss(A, B, eye(4), 0);

    Fd   = exp(-opt.tau * s);                    % camera transport
    Fzoh = exp(-s * opt.Ts / 2);                 % zero-order hold, half-sample
    Fdd  = exp(-s * opt.Ts / 2) / (1 + s/wc);    % backward difference + low-pass

    F = blkdiag(1, 1, Fd, Fd*Fdd);               % [theta, theta_dot, psi, psi_dot]

    L       = K * F * G * Fzoh;                  % implemented, broken at the plant input
    L_ideal = K * G;                             % what the LQR was designed against

    %% --- Margins ---
    % On an frd grid rather than the tf: L carries irrational delay terms,
    % and evaluating them exactly on a grid is honest where a Pade fit
    % would quietly trade away the very phase this function is measuring.
    Lf = frd(L, opt.w);
    [Gm, Pm, Wcg, Wcp] = margin(Lf);

    [Gmi, Pmi, Wcgi, Wcpi] = margin(frd(L_ideal, opt.w));

    out.Gm_dB = 20*log10(Gm);
    out.Pm_deg = Pm;
    out.Wcg = Wcg;
    out.Wcp = Wcp;
    out.ideal = struct('Gm_dB', 20*log10(Gmi), 'Pm_deg', Pmi, 'Wcg', Wcgi, 'Wcp', Wcpi);

    % Extra transport delay the loop still tolerates: at the gain
    % crossover, every added second of delay costs Wcp radians of phase.
    if isfinite(Pm) && Wcp > 0
        out.delay_margin_s = deg2rad(Pm) / Wcp;
    else
        out.delay_margin_s = NaN;
    end

    out.L = L;  out.L_ideal = L_ideal;
    out.K = K;  out.A = A;  out.B = B;
    out.opt = opt;

    %% --- Report ---
    fprintf('\n=== Loop margins: %s, psi_lin = %.1f deg ===\n', ...
        opt.mode, rad2deg(opt.psi_lin));
    fprintf('K  = [%+.4f %+.4f %+.4f %+.4f]\n', K);
    fprintf('Ts = %.1f ms | camera delay = %.1f ms | psi_dot cutoff = %.1f Hz\n', ...
        1e3*opt.Ts, 1e3*opt.tau, opt.fc);
    fprintf('\n%-24s %10s %10s\n', '', 'ideal', 'implemented');
    fprintf('%-24s %10.2f %10.2f\n', 'gain margin [dB]',   out.ideal.Gm_dB,  out.Gm_dB);
    fprintf('%-24s %10.2f %10.2f\n', 'phase margin [deg]', out.ideal.Pm_deg, out.Pm_deg);
    fprintf('%-24s %10.3f %10.3f\n', 'gain crossover [Hz]',  out.ideal.Wcp/(2*pi), out.Wcp/(2*pi));
    fprintf('%-24s %10.3f %10.3f\n', 'phase crossover [Hz]', out.ideal.Wcg/(2*pi), out.Wcg/(2*pi));
    fprintf('\nremaining delay margin: %.1f ms on top of the %.1f ms already modelled\n', ...
        1e3*out.delay_margin_s, 1e3*(opt.tau + opt.Ts));

    if out.Pm_deg < 30
        fprintf(['\nWARNING: %.1f deg of phase margin. Below about 30 deg the step\n' ...
                 'response rings and the loop is one modelling error from instability.\n'], ...
                out.Pm_deg);
    end

    %% --- Plot ---
    if opt.plot
        figure('Name', sprintf('Loop margins -- %s', opt.mode));
        bode(L_ideal, opt.w, 'k--'); hold on;
        bode(L, opt.w, 'b');
        grid on;
        legend('designed (no delay, no filter)', 'implemented', 'Location', 'southwest');
        title(sprintf('%s: Pm %.1f deg -> %.1f deg', ...
            opt.mode, out.ideal.Pm_deg, out.Pm_deg));
    end
end
