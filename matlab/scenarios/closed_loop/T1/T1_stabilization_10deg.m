function run = T1_stabilization_10deg()
    params = ball_hoop_params();
    [~, R_eff_o] = hoop_geometry('rolling_out', params);

    scn.id   = 'T1_stabilization_10deg';
    scn.name = 'T1 - Recovery from 10 deg (outer hoop)';
    scn.mode = 'rolling_out';
    scn.x0   = [R_eff_o; 0; 0; 0; deg2rad(10); deg2rad(0)];

    scn.controller = struct( ...
        'type',    'LQR', ...
        'psi_lin', 0, ...  
        'Q',       diag([1/(10*pi)^2, 1/2^2, 1/0.05^2, 1/0.3^2]), ...
        'R',       0.1);
 
    scn.reference = reference_constant(0);

    run = run_scenario(scn);
end
