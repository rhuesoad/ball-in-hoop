function [hole_start, hole_end] = hole_bounds(params, hoop_theta)
%HOLE_BOUNDS  Angular boundaries of the inner-hoop gap, in the inertial
%             frame, at the current hoop rotation.
%
%   [hole_start, hole_end] = hole_bounds(params, hoop_theta)
%
%   The gap is a cutout in the inner hoop, so its angular position is
%   fixed in the hoop's own (body) frame at params.hole_center_angle, and
%   rotates rigidly with the hoop. psi is measured in the inertial frame
%   (notation contract), so any test of "is the ball over the hole" must
%   compare against the hole's position in the inertial frame, which is
%   hole_center_angle + hoop_theta -- not hole_center_angle alone.
%
%   Inputs
%   ------
%   params     : physical parameters struct (see ball_hoop_params)
%   hoop_theta : current hoop angle [rad]
%
%   Outputs
%   -------
%   hole_start, hole_end : angular boundaries, wrapped to [0, 2*pi) [rad]
%                           Note hole_start > hole_end is possible (the
%                           gap straddles the 0/2*pi wrap); callers must
%                           handle that case, see detect_events.m.

    hole_start = mod(params.hole_center_angle - params.hole_angular_width/2 + hoop_theta, 2*pi);
    hole_end   = mod(params.hole_center_angle + params.hole_angular_width/2 + hoop_theta, 2*pi);
end
