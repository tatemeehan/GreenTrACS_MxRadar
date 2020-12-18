function [ VRMS ] = DryCrimVRMS( rho )
%DryCrimVRMS Returns the RMS Velocity Caused by a SnowPack of Density rho
%  
%
Va = 2.998E8;   % Velocity of Free Space [m/s]
rhoi = 916;       % Pure Ice Density, Ulaby et al.(1986)
Ki = 3.15;      % Real Dielectric Constant of Ice at Microwave Frequency 
                % Ulaby et al. (1986)
                
% Handle Units of Density
if rho<=1
    rho = rho.*1000;
end
                
Phi = 1 - (rho/rhoi);
VRMS = (Va./(Phi.*(1-sqrt(Ki))+sqrt(Ki))).*10^-9;


end

