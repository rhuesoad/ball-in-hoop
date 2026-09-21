function run = OL8_transition_inner_to_outer()
%OL8_TRANSITION_INNER_TO_OUTER  Ball released near the inner hoop with
%                               the hoop itself rotated (theta = pi),
%                               so the gap is positioned to let the ball
%                               transit toward the outer hoop. No
%                               controller: passive dynamics only.
%
%   run = OL8_transition_inner_to_outer()

    params = ball_hoop_params();

    scn.id   = 'OL8_transition_inner_to_outer';
    scn.name = 'Transition from inner to outer hoop';
    scn.mode = 'free_fall';
    scn.x0   = [params.inner_hoop_radius - params.ball_radius - 0.001; 0; pi; 0; pi/3; 0];

    run = run_scenario(scn);
end
