function [s, v, a, j, active, done]= traj4c(t, t_0, s_0, s_end, v_lim, a_lim, j_lim)
% TRAJ4C Calculates a trajectory that is 4 times continuously
% differentiable.
% 
% This function implements a trajectory with sin²-jerk and optional
% constant-acceleration and constant-velocity phases. The duration of
% theses phases is optimized for minimal travel time.
%
% -----
% Author: Matthias Zöbl
% Contact: matthias-zoebl@drei.at
% Version: v1.0
% Created: 17.01.2022, using Matlab R2021b
% -----
% MIT License
% Copyright (c) 2022 Matthias Zöbl
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% -----

%% Preparations
% Compensate for start time and initial position
t = t - t_0;
s_end = s_end - s_0;

% Enforce positive displacement, save sign to invert outputs if required
if s_end < 0
    sign = -1;
    s_end = -s_end;
else
    sign = 1;
end

%% Determine validity of configurations and calculate end time if valid
% Configuration 1: trajectory with constant-velocity phase
if s_end > sqrt(8 * v_lim^3 / min(j_lim, 2 * a_lim^2 / v_lim))
    valid_1 = true;
    t_end_1 = s_end / v_lim + sqrt(8 * v_lim / min(j_lim, 2 * a_lim^2 / v_lim));
else
    valid_1 = false;
    t_end_1 = nan;
end
% Configuration 2: trajectory with constant-acceleration phase
if s_end < 2 * min(v_lim, (sqrt(a_lim^4 + s_end * a_lim * j_lim^2) - a_lim^2) / j_lim)^2 / a_lim
    valid_2 = true;
    t_end_2 = 2 * s_end / min(v_lim, (sqrt(a_lim^4 + s_end * a_lim * j_lim^2) - a_lim^2) / j_lim);
else
    valid_2 = false;
    t_end_2 = nan;
end
% Configuration 3: trajectory with constant-velocity and constant-acceleration phase
if s_end > v_lim^2 / a_lim + 2 * a_lim * v_lim / j_lim && v_lim > 2 * a_lim^2 / j_lim
    valid_3 = true;
    t_end_3 = s_end / v_lim + v_lim / a_lim + 2 * a_lim / j_lim;
else
    valid_3 = false;
    t_end_3 = nan;
end

% Determine fastest valid trajectory configuration
% If none is valid, use configuration 0 (no constant phases)
if valid_1 && (~valid_2 || (t_end_1 < t_end_2)) && (~valid_3 || (t_end_1 < t_end_3))
    conf = 1;
elseif valid_2 && (~valid_1 || (t_end_2 < t_end_1)) && (~valid_3 || (t_end_2 < t_end_3))
    conf = 2;
elseif valid_3 && (~valid_1 || (t_end_3 < t_end_1)) && (~valid_2 || (t_end_3 < t_end_2))
    conf = 3;
else
    conf = 0;
end

% Calculate trajectory parameters acording to selected configuration
if conf == 1
    %v_max = v_lim;
    j_max = min(j_lim, 2 * a_lim^2 / v_lim);
    %a_max = sqrt(j_max * v_lim / 2);
    alpha = sqrt(2 * v_lim / j_max);
    t_a_max = 0;
    t_v_max = s_end / v_lim - sqrt(8 * v_lim / j_max);
elseif conf == 2
    %a_max = a_lim;
    v_max = min(v_lim, (sqrt(a_lim^4 + s_end * a_lim * j_lim^2) - a_lim^2) / j_lim);
    j_max = 2 * a_lim^2 * v_max / (a_lim * s_end - v_max^2);
    alpha = s_end / v_max - v_max / a_lim;
    t_a_max = (2 * v_max^2 - a_lim * s_end) / (a_lim * v_max);
    t_v_max = 0;
elseif conf == 3
    %a_max = a_lim;
    %v_max = v_lim;
    j_max = j_lim;
    alpha = 2 * a_lim / j_lim;
    t_a_max = v_lim / a_lim - 2 * a_lim / j_lim;
    t_v_max = s_end / v_lim - v_lim / a_lim - 2 * a_lim / j_lim;
else
    v_max = min([v_lim, sqrt(a_lim * s_end / 2), (j_lim * s_end^2 / 8)^(1 / 3)]);
    %a_max = 2 * v_max^2 / s_end;
    j_max = 8 * v_max^3 / s_end^2;
    alpha = sqrt(s_end^2 / (4 * v_max^2));
    t_a_max = 0;
    t_v_max = 0;
end

%%  Determine flags
active = 0 <= t && t <= 4 * alpha + 2 * t_a_max + t_v_max;
done = t > 4 * alpha + 2 * t_a_max + t_v_max;

% Calculate values for s, v, a and j
if t < 0
    j = 0;
    a = 0;
    v = 0;
    s = 0;
elseif t <= alpha
    j = j_max * sin(t * pi / alpha)^2;
    a = -j_max * (alpha * sin(2 * pi * t / alpha) - 2 * pi * t) / (4 * pi);
    v = t^2 * j_max / 4 + alpha^2 * j_max * (cos(pi * t / alpha)^2 - 1) / (4 * pi^2);
    s = j_max * (4 * pi^3 * t^3 + 3 * alpha^3 * sin(2 * pi * t / alpha) - 6 * pi * t *alpha^2) / (48 * pi^3);
elseif t <= alpha + t_a_max
    j = 0;
    a = alpha * j_max / 2;
    v = alpha^2 * j_max / 4 + alpha * j_max * (t - alpha) / 2;
    s = (2 * pi^2 - 3) / (24 * pi^2) * alpha^3 * j_max + t * alpha * j_max * (t - alpha) / 4;
elseif t <= 2 * alpha + t_a_max
    j = -j_max * sin((t - alpha - t_a_max) * pi / alpha)^2;
    a = alpha * j_max / 2 + j_max * (2 * pi * (alpha - t + t_a_max) + alpha * sin((2 * pi * (t - t_a_max)) / alpha)) / (4 * pi);
    v = ((1 / 4 - (2 * cos(pi * (t - t_a_max) / alpha)^2 + 6 * pi^2 - 2) / (8 * pi^2)) * alpha^2 + (t - t_a_max / 2) * alpha + t * t_a_max / 2 - t_a_max^2 / 4 - t^2 / 4) * j_max;
    s = ((-(12 * pi + 3 * sin(2 * pi * (t - t_a_max) / alpha) - 8 * pi^3) / (48 * pi^3)) * alpha^3 + ((12 * pi^3 - 6 * pi) / (48 * pi^3) * t_a_max + (6 * pi * t - 24 * pi^3 * t) / (48 * pi^3)) * alpha^2 + (t_a_max^2 / 4 - t * t_a_max / 2 + t^2 / 2) * alpha + t_a_max^3 / 12 - t * t_a_max^2 / 4 + t^2 * t_a_max / 4 - t^3 / 12) * j_max;
elseif t <= 2 * alpha + t_a_max + t_v_max
    j = 0;
    a = 0;
    v = alpha * j_max * (alpha + t_a_max) / 2;
    s = ((t / 2 - 3 * t_a_max / 4) * alpha^2 - alpha^3 / 2 + (t * t_a_max / 2 - t_a_max^2 / 4) * alpha) * j_max;
elseif t <= 3 * alpha + t_a_max + t_v_max
    j = -j_max * sin((t - 2 * alpha - t_a_max - t_v_max) * pi / alpha)^2;
    a = j_max * (2 * pi * (2 * alpha - t + t_a_max + t_v_max) - alpha * sin((2 * pi * (t_a_max - t + t_v_max)) / alpha)) / (4 * pi);
    v = ((1 / 2 - (cos(pi * (t_a_max - t + t_v_max) / alpha)^2 + 4 * pi^2 - 1) / (4 * pi^2)) * alpha^2 + (t - t_v_max - t_a_max / 2) * alpha + (t / 2 - t_v_max / 2) * t_a_max - t_a_max^2 / 4 + t * t_v_max / 2 - t_v_max^2 / 4 - t^2 / 4) * j_max;
    s = (((3 * sin(2 * pi * (t_a_max - t + t_v_max) / alpha) - 12 * pi + 8 * pi^3) / (48 * pi^3)) * alpha^3 + ((12 * pi^3 - 6 * pi) / (48 * pi^3) * t_a_max + (48 * pi^3 - 6 * pi) / (48 * pi^3) * t_v_max + (6 * pi * t - 24 * pi^3 * t) / (48 * pi^3)) * alpha^2 + (t_a_max^2 / 4 + (t_v_max - t / 2) * t_a_max + t_v_max^2 / 2 - t * t_v_max + t^2 / 2) * alpha + t_a_max^3 / 12 + (t_v_max / 4 - t / 4) * t_a_max^2 + (t_v_max^2 / 4 - t * t_v_max / 2 + t^2 / 4) * t_a_max + t_v_max^3 / 12 - t * t_v_max^2 / 4 + t^2 * t_v_max / 4 - t^3 / 12) * j_max;
elseif t <= 3 * alpha + 2 * t_a_max + t_v_max
    j = 0;
    a = -alpha * j_max / 2;
    v = (7 * alpha^2 / 4 + (t_a_max + t_v_max / 2 - t / 2) * alpha) * j_max;
    s = ((1 / (8 * pi^2) - 25 / 12) * alpha^3 + (7 * t / 4 - 5 * t_v_max / 4 - 2 * t_a_max) * alpha^2 + ((t - t_v_max / 2) * t_a_max - t_a_max^2 / 2 + t * t_v_max / 2 - t_v_max^2 / 4 - t^2 / 4) * alpha) * j_max;
elseif t <= 4 * alpha + 2 * t_a_max + t_v_max
    j = j_max * sin((t - 3 * alpha - 2 * t_a_max - t_v_max) * pi / alpha)^2;
    a = -alpha * j_max / 2 - j_max * (2 * pi * (3 * alpha - t + 2 * t_a_max + t_v_max) - alpha * sin((2 * pi * (2 * t_a_max - t + t_v_max)) / alpha)) / (4 * pi);
    v = ((1 / 4 - (sin(pi * (2 * t_a_max - t + t_v_max) / alpha)^2 - 15 * pi^2) / (4 * pi^2)) * alpha^2 + (4 * t_a_max + 2 * t_v_max - 2 * t) * alpha + t_a_max^2 + (t_v_max - t) * t_a_max + t_v_max^2 / 4 - t * t_v_max / 2 + t^2 / 4) * j_max;
    s = ((-(3 * sin(2 * pi * (2 * t_a_max - t + t_v_max) / alpha) - 24 * pi + 208 * pi^3) / (48 * pi^3)) * alpha^3 + (-(312 * pi^3 - 12 * pi) / (48 * pi^3) * t_a_max - (168 * pi^3 - 6 * pi) / (48 * pi^3) * t_v_max - (6 * pi * t - 192 * pi^3 * t) / (48 * pi^3)) * alpha^2 + (-7 * t_a_max^2 / 2 + (4 * t - 7 * t_v_max / 2) * t_a_max - t_v_max^2 + 2 * t * t_v_max - t^2) * alpha - 2 * t_a_max^3 / 3 + (t - t_v_max) * t_a_max^2 + (t * t_v_max - t_v_max^2 / 2 - t^2 / 2) * t_a_max - t_v_max^3 / 12 + t * t_v_max^2 / 4 - t^2 * t_v_max / 4 + t^3 / 12) * j_max;
else
    j = 0;
    a = 0;
    v = 0;
    s = (alpha^3 + (3 * t_a_max / 2 + t_v_max / 2) * alpha^2 + (t_a_max^2 / 2 + t_v_max * t_a_max / 2) * alpha) * j_max;
end

% Compensate for sign and initial position
j = j * sign;
a = a * sign;
v = v * sign;
s = s * sign;
s = s + s_0;

end