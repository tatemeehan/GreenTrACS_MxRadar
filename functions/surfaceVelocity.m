function [Vresult] = surfaceVelocity(dataDir,filenameX,filenameY,trhd)
% surfaceVelocity uses NASA MEaSUREs ice sheet surface velocities as a
% first order correction for the MxHL predicted isochrone horizons.
%   The multi-channel gpr array traverses at some angle oblique to the ice
%   sheet motion. This code corrects the ice sheet velocity magnitudes such
%   that they represent the ice displacement velocity on strike with the 
%   GPR array.
%
% This code is designed to work with:
% MEaSUREs Multi-year Greenland Ice Sheet Velocity Mosaic, Version 1
% Joughin, I., B. Smith, I. Howat, and T. Scambos. 2016. 
% Boulder, Colorado USA. NASA National Snow and Ice Data Center Distributed
% Active Archive Center. doi: https://doi.org/10.5067/QUA5Q9SVMSJG
% 
% The user must download then specify the directory (dataDir) containing 
% the vx and vy geotiff files (filenameX, filenameY).
%
% Boise State University: Tate Meehan, NASA ISGC 2019

% Load MEaSUREs Multi-year Greenland Ice Sheet Velocity Mosaic, Version 1. 
[Vx,GeoLoc] = geotiffread(fullfile(dataDir,filenameX));
[Vy,~] = geotiffread(fullfile(dataDir,filenameY));
% Convert no data values to NaN
nanIx = find(Vx<=-2e9);
Vx(nanIx) = NaN;
Vy(nanIx) = NaN;
% Store Array Size
[ny, nx] = size(Vx);
% Extract Grid Axes - NSIDC Sea Ice Polar Stereographic North Projection
dX = GeoLoc.XWorldLimits; dY = GeoLoc.YWorldLimits;
X = linspace(dX(1),dX(2),nx); Y = fliplr(linspace(dY(1),dY(2),ny));
% Read Radar Array GPS from Trace Header
[~, ix] = unique(trhd(10,:));
ix = sort(ix);
Xgps = trhd(10,ix);
Ygps = trhd(11,ix);
heading = trhd(20,ix);

ntrc = length(ix);

% Search Velocity Grid for Radar Array Position
Vix = zeros(ntrc,1);
if isempty(gcp('nocreate'))
    for ii = 1:ntrc
        [~,Xix] = min(abs(X-Xgps(ii)));
        [~,Yix] = min(abs(Y-Ygps(ii)));
        Vix(ii) = sub2ind([ny,nx],Yix,Xix);
    end
else % Use Parallel Loop if worker pool is running
    parfor ii = 1:ntrc
        [~,Xix] = min(abs(X-Xgps(ii)));
        [~,Yix] = min(abs(Y-Ygps(ii)));
        Vix(ii) = sub2ind([ny,nx],Yix,Xix);
    end
end

% Calculate Velocity Magnitude with Euclidian Norm
Vmag = sqrt(Vy(Vix).^2+Vx(Vix).^2);
% atan(X/Y) implies a pi/2 phase rotation for mapping 0 relative to North
Vtheta = atan2d(Vx(Vix),Vy(Vix));
% Convert to 0 - 360
Vtheta = Vtheta + (Vtheta < 0)*360;
% Calculate Resultant Magnitude 
% Direction of Ice Motion Relative to Array Direction
dTheta = heading - Vtheta;
Vresult = Vmag.*cosd(dTheta);

end

