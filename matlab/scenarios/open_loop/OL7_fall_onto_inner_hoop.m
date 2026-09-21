function run = OL7_fall_onto_inner_hoop()
%OL7_FALL_ONTO_INNER_HOOP  Ball released above the inner hoop, falling
%                          onto its outside surface. No controller:
%                          passive dynamics only.
%
%   run = OL7_fall_onto_inner_hoop()

    params = ball_hoop_params();

    scn.id   = 'OL7_fall_onto_inner_hoop';
    scn.name = 'Falling onto inner hoop';
    scn.mode = 'free_fall';
    scn.x0   = [params.inner_hoop_radius - params.ball_radius - 0.001; 0; 0; 0; 11*pi/12; 0];

    run = run_scenario(scn);
end
