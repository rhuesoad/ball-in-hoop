function run = T2_tracking_5deg_1hz()
%T2_TRACKING_5DEG_1HZ  LQR tracks a 5 deg / 1 Hz sinusoid, outer hoop.
%   Matches the bench runs in python/analysis/npz (T2_A05.0deg_f01.00Hz).

    params = ball_hoop_params();
    [~, R_eff_o] = hoop_geometry('rolling_out', params);

    scn.id   = 'T2_tracking_5deg_1hz';
    scn.name = 'T2 - Sine 5 deg at 1 Hz (outer)';
    scn.mode = 'rolling_out';
    scn.x0   = [R_eff_o; 0; 0; 0; 0; 0];

    scn.controller = struct( ...
        'type',    'LQR', ...
        'psi_lin', 0, ...
        'Q',       diag([1e-8, 1e-6, 1/deg2rad(1)^2, 1/0.3^2]), ...   % 1 Hz tuning, as T2_tracking_10deg_1hz
        'R',       5e-3);
    scn.reference = reference_sinusoid(deg2rad(5), 1);

    run = run_scenario(scn);
end
