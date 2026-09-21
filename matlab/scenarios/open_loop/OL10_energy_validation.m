function run = OL10_energy_validation()
%OL10_ENERGY_VALIDATION  Ball rolling on the outer hoop with a high
%                        initial speed (1.1x the minimum apex speed for
%                        a full loop), for visually checking energy
%                        conservation over a long, fast, frictionless-
%                        adjacent passive run. No controller.
%
%   run = OL10_energy_validation()

    params = ball_hoop_params();
    [~, R_eff_o] = hoop_geometry('rolling_out', params);

    scn.id   = 'OL10_energy_validation';
    scn.name = 'Energy-validation';
    scn.mode = 'rolling_out';
    scn.x0   = [R_eff_o; 0; -deg2rad(30); 0; 0; 1.2* sqrt(2.3*params.gravity/R_eff_o)];

    run = run_scenario(scn);
end
