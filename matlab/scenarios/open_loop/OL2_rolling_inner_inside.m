function run = OL2_rolling_inner_inside()
%OL2_ROLLING_INNER_INSIDE  Ball rolling on the inside surface of the
%                          inner hoop, released 15 deg from the bottom.
%                          No controller: passive dynamics only.
%
%   run = OL2_rolling_inner_inside()

    params = ball_hoop_params();
    [~, R_eff_i] = hoop_geometry('rolling_in_inside', params);

    scn.id   = 'OL2_rolling_inner_inside';
    scn.name = 'Rolling on inner hoop (inside)';
    scn.mode = 'rolling_in_inside';
    scn.x0   = [R_eff_i; 0; 0; 0; pi/4; 0];

    run = run_scenario(scn);
end
