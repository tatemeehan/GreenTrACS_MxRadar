function [ V, to, h, G, m ] = VrmsIrls( x,t)
%Dix Performs the Dix Inversion for Depth and Velocity Estimates
%   Utilizes Green's Method of Linearization and OLS Inverse

%% Green Method
% Square x and t
% Ensure Column Vectors
x = x(:); x = x.^2;
t = t(:); t = (t.^2);

% Ordinary Lest Squares Inverse for Layer Velocity
G = zeros(numel(x),2);
G(:,1) = ones(numel(x),1);
G(:,2) = x;
% for i = 1:numel(x)
%     G(i,1) = 1;
%     G(i,2) = x(i);
% end
% Parameters on the Order of one sample
m = irls(G,t,.04,.0055,1,100);
% m = irls(G,t,.04,.000055,1,10000);
% old Paramereters
% m = irls(G,t,.0001,.00001,1,5000);
% m = irls(G,t,.001,.0001,1,5000);

% Compute RMS Velocity & 2-way Intercept Time
Vrms = (sqrt(m(2)))^(-1);
to = sqrt(m(1));

% % Compute the Residual
% r = t - G*m;
% r2  = norm(r);
% % Estimate Data Standard Deviation
% s = r2/sqrt(size(G,1)-size(G,2));
% % Compute the Covariance Matrix
% C = s^2*inv(G'*G);
% % Weights
% W = eye(2)*(1./C);%W = diag(1./C);
% % Weighted Model Gw and data dw
% Gw = G*W; dw = t*W;
%% Dix Inversion for First Interface Reflection
V = Vrms;

% Layer Thickness
h = (to*V)/2;
%% Return Unsquared Original Picked Values for Error Analysis
% G = sqrt(G);
% m = sqrt(m);
% Hols = sqrt(Hols);
end

