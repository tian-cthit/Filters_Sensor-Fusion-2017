function [hx, Hx] = rangeBearingMeasurements(x, s)
%RANGEBEARINGMEASUREMENTS calculates the range and the bearing to the
%position given by the state vector x, from a sensor locateed in s
%
%Input:
%   x           [n x 1] State vector
%   s           [2 x 1] Sensor position
%
%Output:
%   hx          [2 x 1] measurement model evaluated at state x
%   Hx          [2 x n] measurement model Jacobian evaluated at state x
%
% NOTE: the measurement model assumes that in the state vector x, the first
% two states are X-position and Y-position.

% Your code here
hx = [norm(x(1:2)-s);
    atan2(x(2) - s(2), x(1) - s(1))];
n = length(x);
a = x(1) - s(1);
b = x(2) - s(2);
Hx = [a/norm([a;b]), b/norm([a;b]), zeros(1,n-2);
    -b/(a^2 + b^2), a/(a^2 + b^2), zeros(1,n-2)];
end