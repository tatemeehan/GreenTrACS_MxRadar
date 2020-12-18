function [dout] = AGCgain(d,parameters,type);
%GAIN: Gain a group of traces.
%
%  [dout] = gain(d,dt,option1,parameters,option2);
%
%  IN   d(nt,nx):   traces
%       parameters = [agc_gate], length of the agc gate in secs
%       option2 = 0  No normalization
%               = 1  Normalize each trace by amplitude
%               = 2  Normalize each trace by rms value
%
%  OUT  dout(nt,nx): traces after application of gain function
%
%
%  Example:
%
%    d = hyperbolic_events; dout = gain(d,0.004,'agc',0.05,1);
%    wigb([d,dout]);
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: M.D.Sacchi
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%
%     
if nargin == 2
    type = 0;   % Default is no Normalization
end

 [nt,nx] = size(d);

 
% AGC 

  L = parameters(1);
  L = floor(L/2);
  h = triang(2*L+1);

  for k = 1:nx
   aux =  d(:,k);
   e = aux.^2;
   rms = sqrt(conv(e,h,'same')./sum(h));
   epsi = 1.e-10*max(rms);
   op = rms./(rms.^2+epsi);
   dout(:,k) = d(:,k).*op;
   end


 if type==1;                % Normalize by amplitude 

   for k = 1:nx
    aux =  dout(:,k);
    amax = max(abs(aux));
    dout(:,k) = dout(:,k)/amax;
   end

 end


 if type==2;                % Normalize by rms 

   for k = 1:nx
    aux =  dout(:,k);
    amax =  sqrt(sum(aux.^2)/nt);
    dout(:,k) = dout(:,k)/amax;
   end

 end

