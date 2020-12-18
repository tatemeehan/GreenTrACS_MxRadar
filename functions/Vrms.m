function [ V, to, h ] = Vrms( x,t )
%Dix Performs the Dix Inversion for Depth and Velocity Estimates
%   Utilizes Green's Method of Linearization and OLS Inverse

%% Green Method
% Square x and t
% Ensure Column Vectors
x = x(:); x = x.^2;
t = t(:); t = (t.^2);

% Ordinary Lest Squares Inverse for Layer Velocity
G = zeros(numel(x),2);
for i = 1:numel(x)
    G(i,1) = 1;
    G(i,2) = x(i);
end

% Compute OLS Inverse
Hols = inv(G'*G)*G';
% m1 is the Intercept; m2 is Linear Slope
m = Hols*t;

% Compute RMS Velocity & 2-way Intercept Time
Vrms = (sqrt(m(2)))^(-1);
to = sqrt(m(1));

%% Dix Inversion for First Interface Reflection
V = Vrms;

% Layer Thickness
h = (to*V)/2;
%% Return Unsquared Original Picked Values for Error Analysis
% G = sqrt(G);
% m = sqrt(m);
% Hols = sqrt(Hols);
end

