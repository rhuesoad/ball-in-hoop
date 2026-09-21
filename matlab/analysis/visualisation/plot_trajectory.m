function plot_trajectory(params, r, theta, psi, Xe_all, map2plot)
%PLOT_TRAJECTORY  Ball path overlaid on hoop geometry.
%
%   Inputs
%   ------
%   params   : physical parameters from ball_hoop_params()
%   r        : ball radial position history (N x 1)         [m]
%   theta    : hoop angle history (N x 1)                     [rad]
%   psi      : ball angle history (N x 1)                      [rad]
%   Xe_all   : states at mode-transition events (M x 6),
%              [r, r_dot, theta, theta_dot, psi, psi_dot]
%   map2plot : function handle (x_phys, y_phys) -> (x_plot, y_plot),
%              coordinate transform for display orientation
%
%   Outputs: none (creates a figure)

    % No text is drawn on this figure: no title, no axis labels, no legend.
    % The DisplayName properties below are kept anyway, so that restoring
    % the legend is a one-line change rather than a re-labelling job.
    figure('Name', 'Ball Trajectory');
    hold on;
    axis equal;

    theta_circ = linspace(0, 2*pi, 300);

    % - Outer hoop -
    [x_out_outer, y_out_outer] = map2plot(...
        params.outer_hoop_outer_radius * cos(theta_circ + theta(1)), ...
        params.outer_hoop_outer_radius * sin(theta_circ + theta(1)));
    [x_out_inner, y_out_inner] = map2plot(...
        params.outer_hoop_inner_radius * cos(theta_circ + theta(1)), ...
        params.outer_hoop_inner_radius * sin(theta_circ + theta(1)));

    fill([x_out_outer, fliplr(x_out_inner)], ...
         [y_out_outer, fliplr(y_out_inner)], ...
         [0.7 0.7 0.7], 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(x_out_outer, y_out_outer, 'k-', 'LineWidth', 0.7, 'HandleVisibility', 'off');
    plot(x_out_inner, y_out_inner, 'k-', 'LineWidth', 0.7, 'HandleVisibility', 'off');

    % - Inner hoop (with hole) -
    hole_start = params.hole_center_angle - params.hole_angular_width/2;
    hole_end   = params.hole_center_angle + params.hole_angular_width/2;

    mask1 = theta_circ < hole_start;
    mask2 = theta_circ > hole_end;

    [x_in_outer1, y_in_outer1] = map2plot(...
        params.inner_hoop_outer_radius * cos(theta_circ(mask1) + theta(1)), ...
        params.inner_hoop_outer_radius * sin(theta_circ(mask1) + theta(1)));
    [x_in_outer2, y_in_outer2] = map2plot(...
        params.inner_hoop_outer_radius * cos(theta_circ(mask2) + theta(1)), ...
        params.inner_hoop_outer_radius * sin(theta_circ(mask2) + theta(1)));

    [x_in_inner1, y_in_inner1] = map2plot(...
        params.inner_hoop_inner_radius * cos(theta_circ(mask1) + theta(1)), ...
        params.inner_hoop_inner_radius * sin(theta_circ(mask1) + theta(1)));
    [x_in_inner2, y_in_inner2] = map2plot(...
        params.inner_hoop_inner_radius * cos(theta_circ(mask2) + theta(1)), ...
        params.inner_hoop_inner_radius * sin(theta_circ(mask2) + theta(1)));

    plot([x_in_outer1, NaN, x_in_outer2], [y_in_outer1, NaN, y_in_outer2], ...
         'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot([x_in_inner1, NaN, x_in_inner2], [y_in_inner1, NaN, y_in_inner2], ...
         'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    fill([x_in_outer1, fliplr(x_in_inner1)], ...
         [y_in_outer1, fliplr(y_in_inner1)], ...
         [0.7 0.7 0.7], 'EdgeColor', 'none', 'HandleVisibility', 'off');
    fill([x_in_outer2, fliplr(x_in_inner2)], ...
         [y_in_outer2, fliplr(y_in_inner2)], ...
         [0.7 0.7 0.7], 'EdgeColor', 'none', 'HandleVisibility', 'off');

    [x_edge1_out, y_edge1_out] = map2plot(...
        params.inner_hoop_outer_radius * cos(hole_start + theta(1)), ...
        params.inner_hoop_outer_radius * sin(hole_start  + theta(1)));
    [x_edge1_in, y_edge1_in] = map2plot(...
        params.inner_hoop_inner_radius * cos(hole_start + theta(1)), ...
        params.inner_hoop_inner_radius * sin(hole_start + theta(1)));
    plot([x_edge1_out, x_edge1_in], [y_edge1_out, y_edge1_in], ...
         'k-', 'LineWidth', 0.7, 'HandleVisibility', 'off');

    [x_edge2_out, y_edge2_out] = map2plot(...
        params.inner_hoop_outer_radius * cos(hole_end + theta(1)), ...
        params.inner_hoop_outer_radius * sin(hole_end + theta(1)));
    [x_edge2_in, y_edge2_in] = map2plot(...
        params.inner_hoop_inner_radius * cos(hole_end + theta(1)), ...
        params.inner_hoop_inner_radius * sin(hole_end + theta(1)));
    plot([x_edge2_out, x_edge2_in], [y_edge2_out, y_edge2_in], ...
         'k-', 'LineWidth', 0.7, 'HandleVisibility', 'off');

    lim = params.outer_hoop_inner_radius + 0.05;   % 0.05 m: arbitrary display margin

    [x_ball, y_ball] = map2plot(r .* cos(psi), r .* sin(psi));
    plot(x_ball, y_ball, 'b-', 'LineWidth', 2, 'DisplayName', 'Ball path');
    plot(x_ball(1), y_ball(1), 'bo', 'MarkerSize', 8, 'LineWidth', 2, ...
        'DisplayName', 'Starting Position');

    % Direction of travel. A path that doubles back -- which every
    % oscillation does -- draws the same curve whichever way the ball ran
    % along it, so the figure alone cannot say which. The arrowheads are
    % spaced by ARC LENGTH rather than by time, so a fast flight gets as
    % many as a slow roll instead of the samples deciding, and they are
    % drawn with quiver in data coordinates so they stay on the curve when
    % the axes are rescaled (annotation() places arrows in figure
    % coordinates and would not).
    plot_direction_arrows(x_ball, y_ball, 0.03 * 2 * lim);

    if ~isempty(Xe_all)
        re = Xe_all(:, 1);
        psie = Xe_all(:, 5);
        [xe, ye] = map2plot(re .* cos(psie), re .* sin(psie));
        plot(xe, ye, 'ro', 'MarkerSize', 8, 'LineWidth', 2, ...
             'DisplayName', 'State transitions');
    end

    % Frame, ticks, labels and title all removed: the hoops are the scale
    % here, and the figure is read against the geometry rather than against
    % coordinates. The limits are still set, since they fix how much
    % surrounding space is kept once the axes stop being drawn.
    axis([-lim lim -lim lim]);
    axis off;
end


function plot_direction_arrows(x, y, arrow_length)
%PLOT_DIRECTION_ARROWS  Arrowheads along a path, showing which way it was
%                       travelled.
%
%   Heads only, no shafts: the path is already drawn, and a shaft would
%   double the line over part of it, thickening the curve where an arrow
%   sits and reading as a change of the trajectory rather than an
%   annotation on it. Each head is a filled triangle laid down in data
%   coordinates, so it stays on the curve when the axes are rescaled.
%
%   Inputs
%   ------
%   x, y         : path, in the plotted (display) coordinates
%   arrow_length : largest head length allowed, same units as x/y. The head
%                  actually used is the smaller of this and a fraction of
%                  the path's own extent: an oscillation confined to a few
%                  centimetres would otherwise be decorated with heads as
%                  large as the swing itself.
%
%   Outputs: none (draws into the current axes)

    N_ARROWS = 0;

    step = hypot(diff(x), diff(y));
    s    = [0; cumsum(step(:))];
    if s(end) <= 0
        return;   % the ball never moved: nothing to point at
    end

    extent = max(max(x) - min(x), max(y) - min(y));
    arrow_length = min(arrow_length, 0.20 * extent);

    % A path that retraces itself -- again, every oscillation -- brings
    % several of the equally spaced positions back to nearly the same
    % place, at different phases and so pointing opposite ways. Keeping
    % them all produces a blot rather than a direction, so a head is
    % dropped when it would land on one already drawn.
    placed = zeros(0, 2);
    min_separation = 1.5 * arrow_length;

    % Midpoints of N equal arc-length stretches, so no arrow lands on the
    % start marker or on the final sample.
    targets = s(end) * ((1:N_ARROWS) - 0.5) / N_ARROWS;

    for k = 1:numel(targets)
        i = find(s >= targets(k), 1, 'first');

        % Tangent taken over a finite stretch of path, not between two
        % neighbouring samples: near a turning point consecutive samples can
        % be a few microns apart and their difference is then numerical
        % noise, not a direction.
        j = find(s >= s(i) + 0.02 * s(end), 1, 'first');
        if isempty(j), j = numel(s); end
        if j <= i
            i = max(1, i - 1);
            j = min(numel(s), i + 1);
        end

        dx = x(j) - x(i);
        dy = y(j) - y(i);
        d  = hypot(dx, dy);
        if d <= 0
            continue;
        end

        if ~isempty(placed) && ...
           min(hypot(placed(:,1) - x(i), placed(:,2) - y(i))) < min_separation
            continue;
        end
        placed(end+1, :) = [x(i), y(i)]; %#ok<AGROW>

        % Unit tangent, and the normal to it, used as the local frame the
        % triangle is written in.
        ux = dx / d;
        uy = dy / d;

        L = arrow_length;          % tip-to-base length
        W = 0.55 * arrow_length;   % full base width

        % Triangle centred on the path point: tip ahead, base behind.
        tip   = [x(i) + 0.5*L*ux,               y(i) + 0.5*L*uy];
        back  = [x(i) - 0.5*L*ux,               y(i) - 0.5*L*uy];
        left  = [back(1) - 0.5*W*uy,            back(2) + 0.5*W*ux];
        right = [back(1) + 0.5*W*uy,            back(2) - 0.5*W*ux];

        patch('XData', [tip(1), left(1), right(1)], ...
              'YData', [tip(2), left(2), right(2)], ...
              'FaceColor', 'b', 'EdgeColor', 'none', ...
              'HandleVisibility', 'off');
    end
end
