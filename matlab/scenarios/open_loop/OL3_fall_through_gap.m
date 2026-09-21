function run = OL3_fall_through_gap()
%OL3_FALL_THROUGH_GAP  Ball released just inside the inner-hoop gap,
%                       falling through it. No controller: passive
%                       dynamics only.
%
%   run = OL3_fall_through_gap()

    params = ball_hoop_params();

    scn.id   = 'OL3_fall_through_gap';
    scn.name = 'Free fall through the hole';
    scn.mode = 'free_fall';
    % Released 1 mm inside the outer track's ball-centre orbit radius, so
    % the ball starts just clear of contact.
    [~, R_eff_o] = hoop_geometry('rolling_out', params);
    scn.x0   = [R_eff_o - 0.001; 0; 0; 0; ...
                params.hole_center_angle + 0.15; 0];

    run = run_scenario(scn);
end
