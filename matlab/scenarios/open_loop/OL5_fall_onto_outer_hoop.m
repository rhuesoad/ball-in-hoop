function run = OL5_fall_onto_outer_hoop()
%OL5_FALL_ONTO_OUTER_HOOP  Ball released far from the gap, falling back
%                          onto the outer hoop's rolling surface. No
%                          controller: passive dynamics only.
%
%   run = OL5_fall_onto_outer_hoop()

    params = ball_hoop_params();

    scn.id   = 'OL5_fall_onto_outer_hoop';
    scn.name = 'Falling until outer hoop contact';
    scn.mode = 'free_fall';
    % Released 1 mm inside the outer track's ball-centre orbit radius, so
    % the ball starts just clear of contact.
    [~, R_eff_o] = hoop_geometry('rolling_out', params);
    scn.x0   = [R_eff_o - 0.001; 0; 0; 0; 2*pi/3; 0];

    run = run_scenario(scn);
end
