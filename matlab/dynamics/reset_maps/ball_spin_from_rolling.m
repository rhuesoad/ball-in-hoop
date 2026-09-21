function phi_dot = ball_spin_from_rolling(x, params, mode)
%BALL_SPIN_FROM_ROLLING  The ball's own spin, deduced from the rolling
%                        constraint. Valid ONLY while it is in contact.
%
%   phi_dot = ball_spin_from_rolling(x, params, mode)
%
%   Inputs
%   ------
%   x      : full state [r; r_dot; theta; theta_dot; psi; psi_dot]
%   params : physical parameters (see ball_hoop_params)
%   mode   : 'rolling_out' | 'rolling_in_outside' | 'rolling_in_inside'
%
%   Outputs
%   -------
%   phi_dot : ball spin about its own axis                     [rad/s]
%
%   This is docs/MODEL.md eq. (1) differentiated, the same relation
%   contact_forces.m uses to build the tangential force:
%
%       phi_dot = coupling_sign * (R_track/r_roll) * theta_dot
%                 + (R_eff/r_roll) * psi_dot
%
%   Note r_roll, not ball_radius: the ball sits in the vee of two O-rings
%   and rolls on a chord of itself, so it spins faster than a
%   single-contact ball would to cover the same ground (hoop_geometry.m).
%
%   SIGN CONVENTION, and it is not free. On a stationary hoop the
%   relation reduces to phi_dot = (R_eff/r_roll)*psi_dot, i.e.
%   phi_dot = v_t/r_roll with v_t = R_eff*psi_dot the ball-centre
%   tangential speed. That is exactly the convention handle_impact.m's
%   angular-momentum law assumes, so the two agree by construction rather
%   than by a sign chosen to make a case come out right. The check that
%   confirms it: feeding a spin taken from this function back into that
%   law, with the ball already rolling, returns psi_dot unchanged --
%   verified algebraically in handle_impact.m's docstring and numerically
%   in analysis/studies/T4/t4c_impact_check.m.

    [~, R_eff, R_track, coupling_sign, ~, r_roll] = hoop_geometry(mode, params);

    theta_dot = x(4);
    psi_dot   = x(6);

    phi_dot = coupling_sign * (R_track / r_roll) * theta_dot ...
              + (R_eff / r_roll) * psi_dot;
end
