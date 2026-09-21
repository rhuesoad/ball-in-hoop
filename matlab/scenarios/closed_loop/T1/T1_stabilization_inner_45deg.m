function run = T1_stabilization_inner_45deg()
%T1_STABILIZATION_INNER_45DEG  LQR recovery from 45 deg on the inner hoop.
%   Gains are the bench's set B, the one the ten inner-hoop runs in
%   python/analysis/npz were driven with, imposed rather than re-derived so
%   the simulation answers for the same controller the hardware ran.

    params = ball_hoop_params();
    [~, R_eff_i] = hoop_geometry('rolling_in_inside', params);

    scn.id   = 'T1_stabilization_inner_45deg';
    scn.name = 'T1 - Recovery from 45 deg (inner hoop)';
    scn.mode = 'rolling_in_inside';
    scn.x0   = [R_eff_i; 0; 0; 0; deg2rad(45); 0];

    scn.controller = struct( ...
        'type',    'LQR', ...
        'psi_lin', 0, ...
        'K',       [1.0, 4.0, 18.0, 5.5]);   % bench gains, K3/K4 sign-flipped for MATLAB
    scn.reference = reference_constant(0);

    run = run_scenario(scn);
end
