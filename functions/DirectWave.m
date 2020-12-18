function [ V, to] = DirectWave( x,t )
%DirectWave Calculates the Velocity and Intercept Time of the Direct Arrival
% from a CMP Gather using the OLS inverse. 
% Velocity should ~= Free Space Velocity.
%   User inputs Picked Arrival Times, and Offset Array

% Ordinary Lest Squares Inverse for Layer Velocity
t = t(:);
x = x(:);
G = zeros(numel(x),2);
G(:,1) = ones(numel(x),1);
G(:,2) = x;

% Compute OLS Inverse
Hols = inv(G'*G)*G';
% m1 is the Intercept; m2 is Linear Slope
m = Hols*t;

% Compute RMS Velocity & 2-way Intercept Time
V = m(2)^(-1);
to = m(1);

end

