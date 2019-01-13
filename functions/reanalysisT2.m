function [outputArg1,outputArg2] = reanalysisT2(inputArg1,inputArg2)
% reanalysisT2 
%   Detailed explanation goes here
% Greenland Coordinates from Globe
Ix = load('globe2greenland.csv');
lon = ncread(fullfile(dataDir,t2file),'lon');
lat = ncread(fullfile(dataDir,t2file),'lat');
[glon,glat] = meshgrid(lon,lat);
% Convert to -180, 180
glon = mod((glon+180),360)-180;
[x,y] = ll2psn(glat(Ix(:,3)),glon(Ix(:,3)));
% Load t2 data
t2 = ncread(fullfile(dataDir,t2file),'t2');
t2 = permute(t2,[2,1,3]);
% Montly to Average 2m Air Temperaturesfor Reanalysis Period
t2 = mean(t2,3);
% Extract Tempuratures of Greenland only
t2 = t2(Ix(:,3));

% Compute Nearest Coordinate Points


outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

