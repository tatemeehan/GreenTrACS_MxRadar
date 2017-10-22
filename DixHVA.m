function [ Vint, Hint ] = DixHVA( VrmsI, VrmsJ, toI, toJ )
%   DixHVA Computes the Dix Inversion for Interval Velocity and Thickness.
%   This code samples random Jackknifed estimates of Vrms from the Horizon
%   Velocity Analysis (HVA) to create an ensemble of interval velocity and
%   thickness estimates.
%
%           Inputs:     Vrms(I/J) - Vector of RMS Velocities from HVA
%                       to(I/J)   - Vector of Reflected Arrival Zero Times
%                                 - I/J is the upper/lower layer estimates
%
%           Output:     Vint - Vector of Interval Velocities
%                       Hint - Vector of Layer Thicknesses
%
%   Written by: Tate Meehan, Boise State University, GreenTrACS 2016-2017
%
%% Dix Inversion
MC = 50;                       % Monte Carlo Simulations
draws = 25;                     % Integer of Random Samples
nVi = size(VrmsI,1);
nVj = size(VrmsJ,1); % Number of RMS Velocities Simulated HVA
nVix = 1:nVi;  
nVjx = 1:nVj; % Linear Index Array
Vint = zeros(draws,MC);         % Allocation
Hint = zeros(draws,MC);         % Allocation

% Monte Carlo Simulation for Interval Velocites
for ii = 1:MC
    ixMC = randsample(nVix,draws,'true');
    jxMC = randsample(nVjx,draws,'true');
    Vi = VrmsI(ixMC); % Upper Layer RMS velocity
    Vj = VrmsJ(jxMC); % Lower Layer RMS velocity    
    Ti = toI(ixMC);   % Upper Layer intercept time   
    Tj = toJ(jxMC);   % Lower layer intercept time
    % Interval Velocity from Dix's (1955) Solution
        % L1 regularized
        Vint(:,ii) = sqrt(abs(((Vj.^2.*Tj)-(Vi.^2.*Ti))./(Tj-Ti)));
        % Non-regularized
%         Vint(:,ii) = sqrt((((Vj.^2.*Tj)-(Vi.^2.*Ti))./(Tj-Ti)));
    % Layer Thickness
    Hint(:,ii) = Vint(:,ii).*((Tj-Ti))./2;
end
    % Force Columns
    Vint = Vint(:); Hint = Hint(:);
    % Physically Constrained Solution
    iceIx = find(Vint<0.1689 | Vint > 0.25);
    Vint(iceIx) = [];
    Hint(iceIx) = [];
    
    % Remove Outlying Data
    q = quantile(Vint,[0.16,0.84]);
    iqrIx = find(Vint >= q(1) & Vint <= q(2));
    Vint = Vint(iqrIx);
    Hint = Hint(iqrIx);

end

