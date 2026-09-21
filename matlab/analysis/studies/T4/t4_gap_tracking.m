function out = t4_gap_tracking(run_or_path)
%T4_GAP_TRACKING  Hoop reference that keeps the inner-hoop gap centred on
%                  the ball while it crosses the hoop material.
%
%   out = t4_gap_tracking(run)     % a T4 run struct, or a path to one
%
%   The flight is ballistic and free_fall_matrices.m does not couple theta to
%   (r, psi), so the reference is pure geometry: theta = psi -
%   hole_center_angle, held over the crossing window only. The hold law in
%   t4_flying_ball_controller.m is a straight line in t, so that requirement
%   is least-squared onto one; out.residual_deg is what the line cannot follow.
%
%   Only the first inward passage counts: the ball's centre dips below its
%   settling orbit before rising back onto the track, and including the return
%   would fold the landing into the fit.
%
%   Returns theta_hold, theta_dot_hold, residual_deg, window_ms, half_ball_deg.
    if ischar(run_or_path) || isstring(run_or_path)
        S = load(run_or_path);  run = S.run;
    else
        run = run_or_path;
    end

    p  = run.params;
    Rb = p.ball_radius;
    mi = run.mode_intervals;
    if numel(mi) < 2 || ~strcmp(mi(2).mode, 'free_fall')
        error('t4_gap_tracking:no_flight', ...
              'This run has no free_fall interval to read the crossing from.');
    end
    tL = mi(1).t_end;   tC = mi(2).t_end;

    band = [p.inner_hoop_inner_radius - Rb, p.inner_hoop_outer_radius + Rb];
    seg  = find(run.t_s >= tL & run.t_s <= tC);
    rr   = run.x(seg, 1);
    k1 = find(rr <= band(2), 1, 'first');
    k2 = k1 - 1 + find(rr(k1:end) <= band(1), 1, 'first');
    if isempty(k1) || isempty(k2)
        error('t4_gap_tracking:no_crossing', 'The ball never crosses the annulus.');
    end
    w  = seg(k1:k2);

    t_fly = run.t_s(w) - tL;
    psi   = run.x(w, 5);
    rc    = run.x(w, 1);

    % Least-squares line onto theta = psi - hole_center_angle.
    need = psi - p.hole_center_angle;
    A    = [ones(numel(t_fly),1), t_fly];
    c    = A \ need;

    % Wind theta_hold to the turn the hoop is actually on: only theta mod
    % 2*pi is physical, but the hold law compares against the accumulated
    % theta, so a reference several turns away would command a huge slew.
    theta_L = run.x(find(run.t_s >= tL, 1, 'first'), 3);
    c(1) = c(1) + round((theta_L - c(1)) / (2*pi)) * 2*pi;

    out.theta_hold     = c(1);
    out.theta_dot_hold = c(2);
    out.residual_deg   = rad2deg(max(abs(A*(A\need) - need)));
    out.window_ms      = 1e3 * [t_fly(1), t_fly(end)];
    out.half_ball_deg  = rad2deg(asin(min(1, Rb ./ [rc(1), rc(end)])));

    fprintf('\n  crossing        : %.1f to %.1f ms after liftoff (%.1f ms)\n', ...
            out.window_ms(1), out.window_ms(2), diff(out.window_ms));
    fprintf('  r               : %.1f -> %.1f mm\n', 1e3*rc(1), 1e3*rc(end));
    fprintf('  psi             : %.1f -> %.1f deg (arc %.1f deg)\n', ...
            rad2deg(psi(1)), rad2deg(psi(end)), rad2deg(psi(1)-psi(end)));
    fprintf('  ball half-angle : %.1f -> %.1f deg\n', out.half_ball_deg(1), out.half_ball_deg(2));
    fprintf('  fit residual    : %.2f deg\n', out.residual_deg);
    fprintf('\n  paste into the scenario:\n');
    fprintf('    ''theta_hold'',     %.4f, ...\n', out.theta_hold);
    fprintf('    ''theta_dot_hold'', %.4f, ...\n\n', out.theta_dot_hold);
end
