function ref = reference_constant(psi_target)
%REFERENCE_CONSTANT  Constant-setpoint reference (for stabilization scenarios).
%
%   ref = reference_constant(psi_target)
%
%   Inputs
%   ------
%   psi_target : desired ball angle, held for all t [rad]
%
%   Outputs
%   -------
%   ref : struct with fields
%           .psi_ref     : @(t) -> psi_target                [rad]
%           .psi_dot_ref : @(t) -> 0                          [rad/s]

    ref.psi_ref     = @(t) psi_target;
    ref.psi_dot_ref = @(t) 0;
end
