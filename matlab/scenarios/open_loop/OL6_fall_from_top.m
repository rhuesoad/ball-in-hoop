function run = OL6_fall_from_top()
%OL6_FALL_FROM_TOP  Ball released at the gap centre angle itself
%                    (the degenerate edge of the free-fall entry
%                    condition). No controller: passive dynamics only.
%
%   run = OL6_fall_from_top()

    params = ball_hoop_params();

    scn.id   = 'OL6_fall_from_top';
    scn.name = 'Falling from top';
    scn.mode = 'free_fall';
    % Released 1 mm inside the outer track's ball-centre orbit radius, so
    % the ball starts just clear of contact.
    [~, R_eff_o] = hoop_geometry('rolling_out', params);
    scn.x0   = [R_eff_o - 0.001; 0; 0; 0; ...
                params.hole_center_angle + 1e-6; 0];

    run = run_scenario(scn);
end
