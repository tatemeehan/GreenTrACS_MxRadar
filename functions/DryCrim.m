function [ pF,PhiF ] = DryCrim( Vf )
%   DryCrim         Approximates Dry Firn/Sno Density and Porosity 
%                   Using a 2-Phase Mixing Model
%
%   Inputs:         EM Velocity of Firn/Sno,       Vf
%
%   Outputs:        Density of Firn,               pf
%                   Porosity of Firn,              PhiF
%
%   Defined:
%                   EM Velovity of Ice,            Vi
%                   EM Velocity of Air,            Va
%                   Density of Ice,                pi
%   Units are MKS
%
%   Tate Meehan - Boise State University Geophysics - March 7, 2016
%       Adapted From A.P. Annan, S.W. Cosway, and T. Sigurdsson (1994).
%       GPR For Snow Pack Water Content. Proceedings of the Fifth
%       International Conference on Ground Penetrating Radar, June 1994.

Va = 2.998E8;   % Velocity of Free Space [m/s]
Ki = 3.15;      % Real Dielectric Constant of Ice at Microwave Frequency 
                % Ulaby et al. (1986)
Vi = Va/sqrt(Ki);   % Velocity of Ice [m/s]
pi = 917;        % Pure Ice Density, Herron and Langway (1980)

% Condition Units of Velocity [m/ns] to [m/s]
if Vf < 1E8
    Vf = Vf.*(10^9); 
end

Kf = (Va./Vf).^2; % Approximate Dielectric Constant of Firn eqn. E. 82 Ulaby et al. (1986)
           
PhiF = (sqrt(Kf)-sqrt(Ki))./(1-sqrt(Ki)); %   Firn/Sno Porosity (Annan et al. 1994)
% PhiF = ((Vi.*Va)-(Vf.*Va))/(Vf.*(Vi-Va)); 
pF = pi.*(1-PhiF);                        %   Firn/Sno Density
pF = pF/1000;                             %   Firn/Sno Percent Density


end

