function [t2] = reanalysisT2(dataDir,t2file,XY)
% reanalysisT2 imports 2 meter air temperature data from Birkle (2018) and
% locates the nearest average annual temperature estimates to the radar.
% Annual average temperatures are estimated on approximately the same time
% scale as the GreenTrACS ice core record and imaging depth of the radar
% approximately coincides with the time period for reanalysis models.

% Greenland Coordinates from Globe
Ix = load('globe2greenland.txt');
lon = ncread(fullfile(dataDir,t2file),'lon');
lat = ncread(fullfile(dataDir,t2file),'lat');
[glon,glat] = meshgrid(lon,lat);
% Convert to -180, 180
glon = mod((glon+180),360)-180;
[x,y] = ll2psn(glat(Ix(:,3)),glon(Ix(:,3)));

% This code functions smoothly with 2 meter air temperature data from:
% Sean Birkel. 2018. Greenland surface mass balance derived from climate 
% reanalysis models, 1979-2017. Arctic Data Center. doi:10.18739/A2D21RH75.
% Load t2 data
t2 = ncread(fullfile(dataDir,t2file),'t2');
t2 = permute(t2,[2,1,3]);
% Montly to Average 2m Air Temperatures for Reanalysis Period
t2 = mean(t2,3);
% Extract Tempuratures of Greenland only
t2 = t2(Ix(:,3));
% PSN Coordinates Axes for Greenland
Xax = linspace(-652925,879625,numel(unique(Ix(:,1))));
Yax = linspace(-3384425,-632675,numel(unique(Ix(:,2))));
% ReGridding Temperature Field is more accurate(~1 degree C error)
[X,Y] = meshgrid(Xax,Yax);
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
T2 = griddata(double(x),double(y),double(t2),double(X),double(Y),'natural');
clear('t2','lon','lat','Xax','Yax','x','y','glon','glat');
% Compute Nearest Coordinate Points
tmp = zeros(1,size(XY,2));
for ii = 1:size(XY,2)
    [~,t2Ix] = min(sqrt((XY(1,ii)-X(:)).^2+(XY(2,ii)-Y(:)).^2));
    % Extract nearest t2
        tmp(ii) = T2(t2Ix);
end
% Output t2
t2 = tmp;
end

