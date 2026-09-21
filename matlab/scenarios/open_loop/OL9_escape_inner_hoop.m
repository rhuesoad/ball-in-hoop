function run = OL9_escape_inner_hoop()
%OL9_ESCAPE_INNER_HOOP  Ball rolling inside the inner hoop, released
%                       with the hoop rotated (theta = 2*pi/3) so the
%                       gap is positioned for the ball to escape through
%                       it. No controller: passive dynamics only.
%
%   run = OL9_escape_inner_hoop()

    params = ball_hoop_params();
    [~, R_eff_i] = hoop_geometry('rolling_in_inside', params);

    scn.id   = 'OL9_escape_inner_hoop';
    scn.name = 'Escaping from inner hoop';
    scn.mode = 'rolling_in_inside';
    scn.x0   = [R_eff_i; 0; 2*pi/3; 0; pi/3; 0];

    run = run_scenario(scn);
end
