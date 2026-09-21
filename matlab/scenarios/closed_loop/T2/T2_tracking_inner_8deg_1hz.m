function run = T2_tracking_inner_8deg_1hz()
%T2_TRACKING_INNER_8DEG_1HZ  LQR tracks a 8 deg / 1 Hz sinusoid, inner hoop.
%   The bench only ever ran T2 on the outer hoop, so there is no measurement
%   to match here. Gains are the inner-hoop nominal set; the Q/R used on the
%   outer scenarios asks for 1 deg of tracking, which saturates the motor at
%   this radius.

    params = ball_hoop_params();
    [~, R_eff_i] = hoop_geometry('rolling_in_inside', params);

    scn.id   = 'T2_tracking_inner_8deg_1hz';
    scn.name = 'T2 - Sine 8 deg at 1 Hz (inner)';
    scn.mode = 'rolling_in_inside';
    scn.x0   = [R_eff_i; 0; 0; 0; 0; 0];

    scn.controller = struct( ...
        'type',    'LQR', ...
        'psi_lin', 0, ...
        'K',       [1.0066, 6.1964, 87.8026, 19.8888]);
    scn.reference = reference_sinusoid(deg2rad(8), 1);

    run = run_scenario(scn);
end
