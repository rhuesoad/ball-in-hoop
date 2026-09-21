function run = OL1_rolling_outer()
%OL1_ROLLING_OUTER  Ball rolling on the outer hoop, released 15 deg from
%                    the bottom. No controller: passive dynamics only.
%
%   run = OL1_rolling_outer()

    params = ball_hoop_params();
    [~, R_eff_o] = hoop_geometry('rolling_out', params);

    scn.id   = 'OL1_rolling_outer';
    scn.name = 'Rolling on outer hoop';
    scn.mode = 'rolling_out';
    scn.x0   = [R_eff_o; 0; 0; 0; pi/4; 0];

    run = run_scenario(scn);
end
