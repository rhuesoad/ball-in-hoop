function animate_system(t_all, X_all, params, map2plot, varargin)
%ANIMATE_SYSTEM  Real-time playback of a simulated run, resampled to a
%                 fixed animation frame rate.
%
%   Inputs
%   ------
%   t_all    : time vector (N x 1)                                [s]
%   X_all    : state history (N x 6), [r, r_dot, theta, theta_dot, psi, psi_dot]
%              [m, m/s, rad, rad/s, rad, rad/s]
%   params   : physical parameters from ball_hoop_params()
%   map2plot : function handle (x_phys, y_phys) -> (x_plot, y_plot),
%              coordinate transform for display orientation. Empty picks
%              the one the report figures use, model x pointing down.
%
%   Name/value
%   ----------
%   'video' : path to write an MP4 to. Empty (default) only plays back.
%   'fps'   : frame rate, and the resampling rate of the playback. 30 by
%             default, so the file runs at the same speed as the run did.
%
%   Outputs: none (draws an animated figure, and writes the file if asked)
    if nargin < 4 || isempty(map2plot)
        map2plot = @(x, y) deal(y, -x);
    end
    p = inputParser;
    p.addParameter('video', '', @(s) ischar(s) || isstring(s));
    p.addParameter('fps',   30, @isscalar);
    p.parse(varargin{:});
    opt = p.Results;

    writing = ~isempty(char(opt.video));
    if writing
        out_dir = fileparts(char(opt.video));
        if ~isempty(out_dir) && ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end
        vw = VideoWriter(char(opt.video), 'MPEG-4');
        vw.FrameRate = opt.fps;
        vw.Quality   = 100;
        open(vw);
        cleanup = onCleanup(@() close(vw));   % the file stays playable if the loop throws
    end

    % Fixed size and offscreen: every frame must come out the same size, and
    % getframe reads the window, which does not survive a long -batch run.
    fig = figure('Name', 'Ball-in-Hoop Animation', ...
                 'Position', [100 100 640 640], ...
                 'Visible', onoff(~writing), 'Color', 'w');
    hold on;
    axis equal;
    grid on;

    lim = params.outer_hoop_outer_radius + 0.02;   % 0.02 m: arbitrary display margin so the hoop edge isn't clipped
    axis([-lim lim -lim lim]);

    xlabel('y [m]');
    ylabel('x [m]');
    title('Real-Time Ball In Hoop Animation');

    dt_anim = 1 / opt.fps;   % playback is real time: one frame per 1/fps of run time
    t_anim = t_all(1):dt_anim:t_all(end);

    [t_unique, ia] = unique(t_all, 'stable');
    X_unique = X_all(ia, :);
    X_interp = interp1(t_unique, X_unique, t_anim);

    theta_circ = linspace(0, 2*pi, 300);
    hole_start = params.hole_center_angle - params.hole_angular_width/2;
    hole_end   = params.hole_center_angle + params.hole_angular_width/2;

    mask1 = theta_circ < hole_start;
    mask2 = theta_circ > hole_end;

    % --- Initialize plot objects ---
    h_outer_fill = fill(NaN, NaN, [0.7 0.7 0.7], 'EdgeColor', 'none');
    h_outer_ext  = plot(NaN, NaN, 'k-', 'LineWidth', 0.7);
    h_outer_int  = plot(NaN, NaN, 'k-', 'LineWidth', 0.7);

    h_inner_fill1 = fill(NaN, NaN, [0.7 0.7 0.7], 'EdgeColor', 'none');
    h_inner_fill2 = fill(NaN, NaN, [0.7 0.7 0.7], 'EdgeColor', 'none');
    h_inner_ext   = plot(NaN, NaN, 'k-', 'LineWidth', 0.7);
    h_inner_int   = plot(NaN, NaN, 'k-', 'LineWidth', 0.7);
    h_edge1       = plot(NaN, NaN, 'k-', 'LineWidth', 0.7);
    h_edge2       = plot(NaN, NaN, 'k-', 'LineWidth', 0.7);

    theta_ball = linspace(0, 2*pi, 50);
    h_ball = fill(NaN, NaN, 'r', 'EdgeColor', 'k', 'LineWidth', 0.5);

    for i = 1:length(t_anim)
        r_i     = X_interp(i, 1);
        theta_i = X_interp(i, 3);
        psi_i   = X_interp(i, 5);

        % === OUTER HOOP ===
        [x_out_outer, y_out_outer] = map2plot(...
            params.outer_hoop_outer_radius * cos(theta_circ + theta_i), ...
            params.outer_hoop_outer_radius * sin(theta_circ + theta_i));
        [x_out_inner, y_out_inner] = map2plot(...
            params.outer_hoop_inner_radius * cos(theta_circ + theta_i), ...
            params.outer_hoop_inner_radius * sin(theta_circ + theta_i));

        set(h_outer_fill, 'XData', [x_out_outer, fliplr(x_out_inner)], ...
                          'YData', [y_out_outer, fliplr(y_out_inner)]);
        set(h_outer_ext, 'XData', x_out_outer, 'YData', y_out_outer);
        set(h_outer_int, 'XData', x_out_inner, 'YData', y_out_inner);

        % === INNER HOOP ===
        [x_in_outer1, y_in_outer1] = map2plot(...
            params.inner_hoop_outer_radius * cos(theta_circ(mask1) + theta_i), ...
            params.inner_hoop_outer_radius * sin(theta_circ(mask1) + theta_i));
        [x_in_outer2, y_in_outer2] = map2plot(...
            params.inner_hoop_outer_radius * cos(theta_circ(mask2) + theta_i), ...
            params.inner_hoop_outer_radius * sin(theta_circ(mask2) + theta_i));

        [x_in_inner1, y_in_inner1] = map2plot(...
            params.inner_hoop_inner_radius * cos(theta_circ(mask1) + theta_i), ...
            params.inner_hoop_inner_radius * sin(theta_circ(mask1) + theta_i));
        [x_in_inner2, y_in_inner2] = map2plot(...
            params.inner_hoop_inner_radius * cos(theta_circ(mask2) + theta_i), ...
            params.inner_hoop_inner_radius * sin(theta_circ(mask2) + theta_i));

        set(h_inner_fill1, 'XData', [x_in_outer1, fliplr(x_in_inner1)], ...
                           'YData', [y_in_outer1, fliplr(y_in_inner1)]);
        set(h_inner_fill2, 'XData', [x_in_outer2, fliplr(x_in_inner2)], ...
                           'YData', [y_in_outer2, fliplr(y_in_inner2)]);

        set(h_inner_ext, 'XData', [x_in_outer1, NaN, x_in_outer2], ...
                         'YData', [y_in_outer1, NaN, y_in_outer2]);
        set(h_inner_int, 'XData', [x_in_inner1, NaN, x_in_inner2], ...
                         'YData', [y_in_inner1, NaN, y_in_inner2]);

        [x_edge1_out, y_edge1_out] = map2plot(...
            params.inner_hoop_outer_radius * cos(hole_start + theta_i), ...
            params.inner_hoop_outer_radius * sin(hole_start + theta_i));
        [x_edge1_in, y_edge1_in] = map2plot(...
            params.inner_hoop_inner_radius * cos(hole_start + theta_i), ...
            params.inner_hoop_inner_radius * sin(hole_start + theta_i));
        [x_edge2_out, y_edge2_out] = map2plot(...
            params.inner_hoop_outer_radius * cos(hole_end + theta_i), ...
            params.inner_hoop_outer_radius * sin(hole_end + theta_i));
        [x_edge2_in, y_edge2_in] = map2plot(...
            params.inner_hoop_inner_radius * cos(hole_end + theta_i), ...
            params.inner_hoop_inner_radius * sin(hole_end + theta_i));

        set(h_edge1, 'XData', [x_edge1_out, x_edge1_in], ...
                     'YData', [y_edge1_out, y_edge1_in]);
        set(h_edge2, 'XData', [x_edge2_out, x_edge2_in], ...
                     'YData', [y_edge2_out, y_edge2_in]);

        % === BALL ===
        x_center = r_i * cos(psi_i);
        y_center = r_i * sin(psi_i);

        x_ball_circle = x_center + params.ball_radius * cos(theta_ball);
        y_ball_circle = y_center + params.ball_radius * sin(theta_ball);

        [x_b, y_b] = map2plot(x_ball_circle, y_ball_circle);
        set(h_ball, 'XData', x_b, 'YData', y_b);
        drawnow;

        if writing
            writeVideo(vw, print(fig, "-RGBImage", "-r100"));
        end
    end

    if writing
        fprintf('Animation written: %s (%d frames, %g fps)\n', ...
                char(opt.video), numel(t_anim), opt.fps);
    else
        fprintf('Animation complete.\n');
    end
end


function s = onoff(tf)
    if tf, s = 'on'; else, s = 'off'; end
end
