function [vH, phiVect,phiMet,u10,v10,U10,V10, stdVect, stdMet] = reanalysisV10U10(dataDir,u10file,v10file,XY)
% reanalysisT2 imports 10 meter v and u wind data from Birkle (2018) and
% locates the nearest average annual velocity estimates to the radar.
% Annual average wind velocities are estimated on approximately the time
% scale of the GreenTrACS ice core record and imaging depth of the radar
% approximately coincides with the time period for reanalysis models.

% Greenland Coordinates from Globe
Ix = load('globe2greenland.txt');
lon = ncread(fullfile(dataDir,u10file),'lon');
lat = ncread(fullfile(dataDir,u10file),'lat');
[glon,glat] = meshgrid(lon,lat);
indx = sub2ind(size(glon),Ix(:,2),Ix(:,1));
% Convert to -180, 180
glon = mod((glon+180),360)-180;
[x,y] = ll2psn(glat(indx),glon(indx));

% This code functions smoothly with 10 meter wind velocity data from:
% Sean Birkel. 2018. Greenland surface mass balance derived from climate 
% reanalysis models, 1979-2017. Arctic Data Center. doi:10.18739/A2D21RH75.
% Load Zonal Wind data
u10 = ncread(fullfile(dataDir,u10file),'u10');
u10 = permute(u10,[2,1,3]);
% Monthly to Average 10 m Zonal Winds for Reanalysis Period
stdU10 = std(u10,[],3);
u10 = mean(u10,3);
% Extract Zonal Winds of Greenland only
u10 = u10(indx);
stdU10 = stdU10(indx);
% PSN Coordinates Axes for Greenland
Xax = linspace(-652925,879625,numel(unique(Ix(:,1))));
Yax = linspace(-3384425,-632675,numel(unique(Ix(:,2))));
% ReGridding Wind Field is more accurate
[X,Y] = meshgrid(Xax,Yax);
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
U10 = griddata(double(x),double(y),double(u10),double(X),double(Y),'natural');
U10std = griddata(double(x),double(y),double(stdU10),double(X),double(Y),'natural');

% Load Zonal Wind data
v10 = ncread(fullfile(dataDir,v10file),'v10');
v10 = permute(v10,[2,1,3]);
% Monthly to Average 10 m Zonal Winds for Reanalysis Period
stdV10 = std(v10,[],3);
v10 = mean(v10,3);
% Extract Zonal Winds of Greenland only
v10 = v10(indx);
stdV10 = stdV10(indx);
% ReGridding Wind Field is more accurate
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
V10 = griddata(double(x),double(y),double(v10),double(X),double(Y),'natural');
V10std = griddata(double(x),double(y),double(stdV10),double(X),double(Y),'natural');

clear('u10','v10','stdV10','stdU10','lon','lat','Xax','Yax','x','y','glon','glat');
% Compute Nearest Coordinate Points
tmpU = zeros(1,size(XY,2));
tmpV = tmpU;
for ii = 1:size(XY,2)
    [~,wIx] = min(sqrt((XY(1,ii)-X(:)).^2+(XY(2,ii)-Y(:)).^2));
    % Extract nearest u10
        tmpU(ii) = U10(wIx);
        tmpV(ii) = V10(wIx);
        tmpUstd(ii) = U10std(wIx);
        tmpVstd(ii) = V10std(wIx);
end
% Output u10 and v10 winds
u10 = tmpU;
u10std = tmpUstd;
v10 = tmpV;
v10std = tmpVstd;
% Compute Resultant Horizontal Wind Velocity 
vH = sqrt(u10.^2+v10.^2);
% Compute Vecotral Wind Direction
phiVect = atan2d(u10,v10);
% Compute Meterorlogical Wind Direction
phiMet = phiVect + 180;
% Convert to 0 - 360 Azimuths
phiVect = phiVect + (phiVect < 0).*360;
phiMet = phiMet + (phiMet < 0).*360;
% Monte Carlo Propoagtion for Uncertainty in Direction
nMC = 1000; pMC = 0.1;
for ii = 1:nMC
    ix = datasample(1:length(u10),round(length(u10).*pMC));
    vHmc(:,ii) = sqrt(u10(ix).^2+v10(ix).^2);
    % Compute Vecotral Wind Direction
phiVectMC(:,ii) = atan2d(u10(ix),v10(ix));
% Compute Meterorlogical Wind Direction
phiMetMC(:,ii) = phiVect(:,ii) + 180;
% Convert to 0 - 360 Azimuths
phiVectMC(:,ii) = phiVect(:,ii) + (phiVect(:,ii) < 0).*360;
phiMetMC(:,ii) = phiMet(:,ii) + (phiMet(:,ii) < 0).*360;
end

end

