function [ V, to] = DirectWaveIrls( x,t )
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

% Compute IRLS Inverse
% Parameters on the Order of One Sample
m = irls(G,t,.2,.02,1,100);
% Old Parameters
% m = irls(G,t,.1,1,1,500);
% Compute RMS Velocity & 2-way Intercept Time
V = m(2)^(-1);
to = m(1);

end

