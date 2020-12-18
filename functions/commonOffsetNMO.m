function [RadarNMO,xStack,tStack,vStack,M] = ...
    commonOffsetNMO(Radar,dt,f0,offset,TraverseX,StackingTime,StackingVelocity,max_stretch,isInterpolate)
%   NMO: A program for NMO correction.
%
% [RadarNMO,xStack,tStack,vStack,M] = nmo(d,dt,h,tnno,vnmo,max_streatch);
%
%  IN   Radar(nt,nx):   data (common-offset gather)
%       dt:             sampling interval in secs
%       h(nh):          vector of offsets in meters
%       TraverseX:      2D Grid of X Position, interceptTime, and Velocity
%       StackingTime:   
%       StackingVelocity:      
%       max_stretch:    maximum stetch allowed in %
%       isInterpolate:  flag for interpolation
%
%  OUT  RadarNMO:  data after NMO correction
%       xStack:    Gridded X position after interpolation
%       tStack:    Gridded intercept Time after interpolation
%       vStack:    Gridded Stacking Velocities after iterpolation
%       M(nt):     number of x-samples that survived muting
%                  at each time position  
%
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
% Modified By: Tate Meehan, Boise State University, GreenTrACS 2017
% Rather than Correcting MidPoint Gathers from a 1D Velocity Model,
% this Program performs the NMO correction for pre-stacked common-offset
% gathers using a 2-D Stacking Velocity Model.



% Post NMO Correction Filtering
    % Convert Filter Units
    f0Hz = f0 * 1e6;        % [Hz]
    dtSec = dt * 1e-9;      % [s]
    
    % Median Filter Parameters
    % Rank of Median Subtraction Filter
    % Nominal Frequency Pass
    MedFiltR = 2.*(ceil(1/((f0Hz)*dtSec))-1)+1;
    % Low Pass
    %          MedFiltR = 2.*(ceil(1/((f0Hz/2.5)*dtSec))-1)+1;
    % High Pass
    %           MedFiltR = 2.*(ceil(1/((f0Hz*2.5)*dtSec))-1)+1;

    % BandPass Filter Parameters
    fMin = f0Hz/2; % [Hz]
    fMax = f0Hz*2; % [Hz]
    
    % Stack Filter Parameters
    StakFiltR = 10; % Filter Rank

% Intepolate t0,v pairs 

  [nt,nx] = size(Radar);
  % Check Flag Condition
  if isInterpolate ~= 1
      isInterpolate = 0;
  end
  
  if isInterpolate
      % Interpolate Model if is intial Loop
      xStack = ones(nt,1)*TraverseX(1,:);
      tStack = ((0:1:nt-1)*dt)'*ones(1,nx);
      vStack = griddata(TraverseX,StackingTime,StackingVelocity,xStack,tStack,'linear');
      vStack(1,:) = vStack(2,:);
  else
      % Else Interpolation has been Gridded on Previous Loop
      xStack = TraverseX; tStack = StackingTime; vStack = StackingVelocity;
  end
  % Allocate Data Output
  RadarNMO = zeros(size(Radar));
  M = zeros(nt,1);
  
  % Test Units of Allowable Stretch Percentage 
  if max_stretch < 1
      max_stretch = max_stretch*100;
  end
  
  % Loop Through Traces in Spatial Domain
  for ix = 1:nx
      % Loop Through Intercept Times
      for it = 1:nt
          % Calulate the Reflection Travel Time for a fixed Offset
          arg = ( tStack(it,ix)^2 + (offset/vStack(it,ix)).^2 );          
          time = sqrt(arg);
          % Apply Stretch Mute
          stretch = (time-tStack(it,ix))/(tStack(it,ix)+1e-10);
          if stretch<max_stretch/100
              
              M(it)= M(it) + 1;
              
              % Grab Sample Index along the NMO Trajectory
              its = time/dt+1;
              % Round Sample UP and DOWN: a is a weight for interpolation
              it1 = floor(time/dt+1);
              it2 = it1+1;
              a = its-it1;
              
              if it2 <= nt
                  % Remap the Data Along the NMO Trajectory to Zero-Offset
                  RadarNMO(it,ix) = (1-a)*Radar(it1,ix)+a*Radar(it2,ix); 
              end
              
          end
      end
%         % Trace Median Subtraction - Noise Supression
%   RadarNMO(:,ix) = medfilt1( RadarNMO(:,ix), MedFiltR, [], 2,'omitnan','truncate' );
%   
%   % BandPass Filter and Mean Subtraction
%   RadarNMO(:,ix) = filter_dylan( RadarNMO(:,ix), dtSec, fMin, fMax, 8 );
  end
  
  % Trace Median Subtraction - Noise Supression
  RadarNMO = medfilt1( RadarNMO, MedFiltR, [], 2,'omitnan','truncate' );
  
  % BandPass Filter and Mean Subtraction
  RadarNMO = filter_dylan( RadarNMO, dtSec, fMin, fMax, 8 );
  
  % Trace Normalize
  tmp = 2.*((RadarNMO(:)-min(RadarNMO(:)))./(max(RadarNMO(:))-min(RadarNMO(:))))-1;
    dq = quantile(tmp(:),[0.425,.5,.575]);
  tmp = tmp - dq(2);
  % reMap Min and Max values to -1 and 1
  maxIx = find(tmp>1);
  minIx = find(tmp<-1);
  if isempty(maxIx) 
      [~, maxIx] = max(tmp);
      tmp(maxIx) = 1;
  else
      tmp(maxIx) = 1;
  end 
  
  if isempty(minIx)
      [~, minIx] = min(tmp);
      tmp(minIx) = -1;
  else
      tmp(minIx) = -1;
  end
  % RadarNMO is Normalized and Reduced to Zero Median
  RadarNMO = reshape(tmp,size(RadarNMO));
  clear('tmp')
  
  % Trace Stacking
  RadarNMO = StakR(RadarNMO,StakFiltR);

 return
    
