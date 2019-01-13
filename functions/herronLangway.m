function [ pz, tz ] = herronLangway( Z, T, po, A, Yr)
% herronLangway estimates the firn density and age at depth. Calculations 
% are published within Herron and Langway 1980, Firn Densification: 
% An Emperical Model. A vector of depths for each estimate must be
% supplied. Mean annual Temperature, Surface Density, and Mean Annual
% Accumulation are optional suppliments. Rough estimates of -20c, 300
% kgm-3, and 0.5 m.w.e. are the default inputs if not specified.
%
%   Inputs:     Z  -  Depths of the estimated densities [m]
%               T  -  The Mean Annual Air tempurature   [K]
%               po -  The Surface Density               [mg/cm3]
%               A  -  The Mean Annual Snow Accumualtion [m.w.e.]
%               Yr -  The Datum Year of the Snow Pack, Drilling, or Survey
%
%   Output:     pz -  The Herron-Langway Estimate of Density at Depth [mg/cm3]
%               tz -  The Predicted Age-Depth Scale
%
%   Written by: Tate Meehan, Boise State University, GreenTrACS 2017
%   Reference: Herron and Langway 1980, Firn Densification: An Emperical Model

if nargin<2
    T = -20;
    po = 0.3;
    A = 0.5;
    Yr = 0;
end
if nargin<3
    po = 0.3;
    A = 0.5;
    Yr = 0;
end
if nargin < 4
    A = 0.5;
    Yr = 0;
end
if nargin < 5
    Yr = 0;
end

% Convert Celcius to Kelvin
if T < 100
    T = T + 273.15;
end
% Check Units of Density
if po > 1
    po = po./1000;
end

% Force Columns
Z = Z(:);
% Boltzman Ideal Gas Constant
R = 8.314;
% Empirical Constant for p < 0.55
ko = 11.*exp(-10160./(R.*T));
% Empirical Constant for p > 0.55
k1 = 575.*exp(-21400./(R.*T));
% Pure Ice Density
pi = 0.917;
% Density-Depth Term p < 0.55
Zo = exp((pi.*ko.*Z)+log(po./(pi-po)));
% Density Depth Relationship p < 0.55
pz1 = (pi.*Zo)./(1+Zo);
% Age Depth Relationship p < 0.55
tz1 = (1./(ko.*A)).*(log((pi-po)./(pi-pz1)));
% Depth at Critical Density
z55Ix = find(pz1<=0.55,1,'last');
z55 = Z(z55Ix);
t55 = tz1(z55Ix);
% Density-Depth Term p > 0.55
Z1 = exp((pi.*k1.*(Z-z55))./sqrt(A) + log(0.55./(pi-0.55)));
% Density Depth Relationship p > 0.55
pz2 = (pi.*Z1)./(1+Z1);
tz2 = (1./(k1.*sqrt(A))).*(log((pi-0.55)./(pi-pz2))) + t55;
z55Ix2 = find(pz2>0.55,1);
% Depth Density Profile
pz = [pz1(1:z55Ix);pz2(z55Ix2:end)];
% Depth Age Profile
tz = [tz1(1:z55Ix);tz2(z55Ix2:end)];
% Convert Age to Years
tz = abs(tz - Yr);
end