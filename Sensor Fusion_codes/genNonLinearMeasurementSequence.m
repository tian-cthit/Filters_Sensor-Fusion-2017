function Y = genNonLinearMeasurementSequence(X, S, h, R)
%GENNONLINEARMEASUREMENTSEQUENCE generates ovservations of the states 
% sequence X using a non-linear measurement model.
%
%Input:
%   X           [n x N+1] State vector sequence
%   S           [n x N] Sensor position vector sequence
%   h           Measurement model function handle
%   R           [m x m] Measurement noise covariance
%
%Output:
%   Y           [m x N] Measurement sequence
%

% Your code here
N = length(X)-1;
Y = zeros(length(R),N);
for i = 1:N
    Y(:,i) = h(X(:,i+1),S(:,i)) + mvnrnd(zeros(length(R),1),R)';
end

end