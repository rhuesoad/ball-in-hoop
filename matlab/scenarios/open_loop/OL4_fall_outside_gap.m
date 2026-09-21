function run = OL4_fall_outside_gap()
%OL4_FALL_OUTSIDE_GAP  Ball released outside the inner-hoop gap's
%                       angular width, so it falls onto the inner hoop
%                       rather than through it. No controller: passive
%                       dynamics only.
%
%   run = OL4_fall_outside_gap()

    params = ball_hoop_params();

    scn.id   = 'OL4_fall_outside_gap';
    scn.name = 'Falling outside the hole';
    scn.mode = 'free_fall';
    % Released 1 mm inside the outer track's ball-centre orbit radius, so
    % the ball starts just clear of contact.
    [~, R_eff_o] = hoop_geometry('rolling_out', params);
    scn.x0   = [R_eff_o - 0.001; 0; 0; 0; ...
                params.hole_center_angle + 0.6; 0];

    run = run_scenario(scn);
end
