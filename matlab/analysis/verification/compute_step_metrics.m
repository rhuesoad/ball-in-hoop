function [t_rise, t_peak, M_p, t_settle] = compute_step_metrics(t, y, y_target)
%COMPUTE_STEP_METRICS  Rise time, peak time, overshoot, and settling time
%                      for a step response.
%
%   [t_rise, t_peak, M_p, t_settle] = compute_step_metrics(t, y, y_target)
%
%   Inputs
%   ------
%   t        : time vector [s]  (N x 1)
%   y        : output signal (N x 1)
%   y_target : target (final) value [rad or any consistent unit]
%
%   Outputs
%   -------
%   t_rise  : 10%-to-90% rise time [s].  NaN if response never reaches 90%.
%   t_peak  : time at which y reaches its absolute maximum [s]
%   M_p     : percent overshoot relative to step magnitude [%].
%             Negative means undershoot.  NaN if undefined.
%   t_settle: last time at which |y - y_target| > 2% * |step_mag| [s].
%             Returns t(1) if the response is within the band from the start.
%
%   Definitions follow Ogata (2010, Modern Control Engineering, sec. 5-3).

t = t(:);
y = y(:);
if isempty(t) || isempty(y) || numel(t) ~= numel(y)
    t_rise = NaN;  t_peak = NaN;  M_p = NaN;  t_settle = NaN;
    return;
end

y0       = y(1);
step_mag = y_target - y0;
denom    = max(abs(step_mag), eps);     % guards against zero-step division

%% --- Rise time (10% to 90% of step) ---
low_threshold  = y0 + 0.1 * step_mag;
high_threshold = y0 + 0.9 * step_mag;

if step_mag >= 0
    idx_low  = find(y >= low_threshold,  1, 'first');
    idx_high = find(y >= high_threshold, 1, 'first');
else
    idx_low  = find(y <= low_threshold,  1, 'first');
    idx_high = find(y <= high_threshold, 1, 'first');
end

if ~isempty(idx_low) && ~isempty(idx_high) && idx_high > idx_low
    % Linear interpolation for sub-sample accuracy.
    if idx_low > 1
        t_low = interp1(y(idx_low-1:idx_low), t(idx_low-1:idx_low), ...
                        low_threshold, 'linear', 'extrap');
    else
        t_low = t(1);
    end
    if idx_high > 1
        t_high = interp1(y(idx_high-1:idx_high), t(idx_high-1:idx_high), ...
                         high_threshold, 'linear', 'extrap');
    else
        t_high = t(idx_high);
    end
    t_rise = t_high - t_low;
else
    t_rise = NaN;
end

%% --- Peak time and overshoot ---
[y_peak, idx_peak] = max(y);
t_peak = t(idx_peak);
M_p    = (y_peak - y_target) / denom * 100;    % positive = overshoot [%]

%% --- Settling time (2% band) ---
band    = 0.02 * denom;
outside = abs(y - y_target) > band;
if any(outside)
    t_settle = t(find(outside, 1, 'last'));
else
    t_settle = t(1);    % within band from the very first sample
end
end
