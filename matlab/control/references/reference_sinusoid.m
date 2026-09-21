function ref = reference_sinusoid(amplitude, frequency_hz)
%REFERENCE_SINUSOID  Sinusoidal oscillation around zero.
%
%   ref = reference_sinusoid(amplitude, frequency_hz)
%
%   Inputs
%   ------
%   amplitude    : oscillation amplitude [rad]
%   frequency_hz : oscillation frequency [Hz]
%
%   Outputs
%   -------
%   ref : struct with fields
%           .psi_ref     : @(t) -> amplitude*sin(2*pi*f*t)              [rad]
%           .psi_dot_ref : @(t) -> 2*pi*f*amplitude*cos(2*pi*f*t)        [rad/s]
%
%   No .psi_ddot field: see archive/reference_quintic.m and docs/MODEL.md sec. 6.3
%   for why T2's scenarios use the static feedforward.

    ref.psi_ref     = @(t)              amplitude              * sin(2*pi*frequency_hz*t);
    ref.psi_dot_ref = @(t) 2*pi*frequency_hz * amplitude        * cos(2*pi*frequency_hz*t);
end
