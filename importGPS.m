%% importGPS - Loads GreenTrACS GPS Data and Extracs Positions for GPR files.
% This program opens files of three formats. GeoHx and NetR8 are .csv formats
% for 2016 traverse. RTKgps reads .txt format from RTKLIB post-processing.
% Positions are converted from Lat, Lon to PSN. Elevation is corrected from
% WGS84 to Geoid via EGM08. A static shift is applied to correct for antenna 
% height above GrIS surface. Slope, Speed, and Heading are calculated. 
% Data is stored in .mat structure "GeoLocation" to be read into MxRadar.m

%  Boise State University: Tate Meehan GreenTrACS 2017

% clear; close all; clc;
addpath './functions'
%% Load GPS.csv
workingDir = pwd;
% dataDir = '/home/tatemeehan/GreenTracs2017/GPS/Core15SpurW061317';
dataDir = 'D:\GreenTrACS\2017\GPS\Processed\06\17';
writeDir = '/home/tatemeehan/GreenTracs2017/GPS/Core15SpurW061317';
isWrite = 0;
isNetR8 = 0;    % 5 column Format
isGeoHx = 0;    % Core 2 & Core 5 Spirals
isRTKgps = 1;   % RTKgps Format
% Using uiImport
% cd into Data Directory
cd(dataDir)

if isRTKgps
    [filename, pathname] = uigetfile('*.txt');
    opts = detectImportOptions(filename);
else
    [filename, pathname] = uigetfile('*.csv');
end
cd(workingDir)
fileIDs = fullfile(pathname, filename);
if isNetR8
formatSpec = '%f%f%f%{yyyy/MM/dd}D%{HH:mm:ss.SSS}D';
LLZDT = readtable(fileIDs,'Delimiter',',','Format',formatSpec);
end
if isGeoHx
formatSpec = '%f%f%f%{yyyy/MM/dd}D%{HH:mm:ss}D';
LLZDT = readtable(fileIDs,'Delimiter',',','Format',formatSpec);
end
if isRTKgps
    opts.DataLine = 11;
    fields = readtable(fileIDs,opts);
    % Parse to Five Columns
    fields(:,6:end) = [];
    % Sort Table(lon, lat, elv, date, time)
    LLZDT = fields(:,[4,3,5,1,2]);
end
LLZDT = table2cell(LLZDT);

% Place LLZ in Array
Date = [LLZDT{:,4}];
Date = Date(:);
Time = [LLZDT{:,5}];
Time = Time(:);
Lat = [LLZDT{:,2}];
Lat = Lat(:);
Lon = [LLZDT{:,1}];
Lon = Lon(:);
Zwgs84 = [LLZDT{:,3}];
Zwgs84 = Zwgs84(:);

% Remove StaticGPS Traces
% Include Useage of GPS Edges.
win = 11;
[GPSix,GPSixEdges] = removeStaticGPS(Lon,Lat,win);
% GPSix = 1:length(LLZDT);

% Place LLZDT in Array
Date = Date(GPSix);
Date = Date(:);
Time = Time(GPSix);
Time = Time(:);
Lat = Lat(GPSix);
Lat = Lat(:);
Lon = Lon(GPSix);
Lon = Lon(:);
Zwgs84 = Zwgs84(GPSix);
Zwgs84 = Zwgs84(:);
% Correct Elevations to Surface Height from Antenna Base Plane ~192 cm
% (Personal Communication with Bob Hawley, 2019)
Zwgs84 = Zwgs84-1.92; 

% Smooth Lon, Lat, Elevation
R = 25; % Smoothing Window Rank
X = 1:length(Zwgs84); % Estimtes to Smooth

% % Smooth Coordinate Estimate
% Lon = nonParametricSmooth(X,Lon,X,R);
% Lat = nonParametricSmooth(X,Lat,X,R);

% Smooth Elevation Estimate
if isGeoHx || isRTKgps
    kps = 15./3600;
    km = 0.5;
    R = 2.*(round(km./kps./2))+1;
% R = 121; % ~0.5 km Smoothing at 15 kph
end
Zwgs84 = nonParametricSmooth(X,Zwgs84,X,R);

% Geoid Height: EGM2008
N = geoidheight(Lat,Lon,'EGM2008');
N = N(:);

% Correct Elevation to MSL
Zegm08 = Zwgs84-N;

% Convert Lat Lon Elv to NSIDC Sea Ice Polar Stereographic North Projection
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. (2017). Arctic Mapping Tools
[X, Y] = ll2psn(Lat,Lon);
dX = [0;diff(X)];
dY = [0;diff(Y)];
dZ = diff(Zegm08);
dZ = [dZ(1);dZ];

% Compute Distance
dS = sqrt(dX.^2+dY.^2+dZ.^2);
Distance = cumsum(dS); %[m]
DistanceKm = Distance./1000; % [km]

% Compute Velocity
dT = [1;diff(Time)];
dT = seconds(dT);
Velocity = dS./dT; % [m/s]
Velocity = medfilt1(Velocity,R);

% Heading
dH = [0,0;diff([X,Y])];
heading = mod(atan2(dH(:,1), dH(:,2))*180/pi, 360);
heading = medfilt1(heading,R);
tailing = mod(heading+180,360);

% Slope
Slope = atand(dZ./dS);
Slope = medfilt1(Slope,R);

% Limit Lat and Lon for GeoShow
LatLim = [min(Lat) max(Lat)];
LonLim = [min(Lon) max(Lon)];
[LatLim, LonLim] = bufgeoquad(LatLim, LonLim, 0.005, 0.005);

%% Create GPS Structure
% Select Valid Traces for Export
Selection = (1:length(Lat));


    field1 = 'X';               value1 = X(Selection);
    field2 = 'Y';               value2 = Y(Selection);
    field3 = 'Z_EGM08';         value3 = Zegm08(Selection);
%     field4 = 'Zone';            value4 = {Zone(Selection)};
    field5 = 'Lon';             value5 = Lon(Selection);
    field6 = 'Lat';             value6 = Lat(Selection);
    field7 = 'Z_WGS84';         value7 = Zwgs84(Selection);
    field8 = 'N_EGM08';         value8 = N(Selection);
    field9 = 'Date';            value9 = {LLZDT(Selection,4)};
    field10 = 'Time';           value10 = {LLZDT(Selection,5)};
%     field11 = 'GPSix';          value11 = GPSix(Selection);
    field12 = 'GPSixEdges';     value12 = GPSixEdges;
    field13 = 'Slope';          value13 = Slope;
    field14 = 'Speed';          value14 = Velocity;
    field15 = 'Heading';        value15 = heading;
    field16 = 'dX';             value16 = dX;
    field17 = 'dY';             value17 = dY;
    field18 = 'dZ';             value18 = dZ;
    
    GeoLocation = struct(field1,value1,field2,value2,field3,value3,...%field4,value4,...
        field16,value16,field17,value17,field18,value18,...
        field5,value5,field6,value6,field7,value7,field8,value8,...
        field13,value13,field14,value14,field15,value15,...
        field9,value9,field10,value10,...
        field12,value12);%,field11,value11);


% Save GPS Data
if isWrite
cd(writeDir)
writeName = 'MxRadarGPSCore15SpurW061317.mat';
save(writeName,'GeoLocation','-v7.3')
cd(workingDir)
end

%% Create Maps

figure(1);
ax = worldmap('Greenland');
greenland = shaperead('landareas', 'UseGeoCoords', true,...
  'Selector',{@(name) strcmp(name,'Greenland'), 'Name'});
patchm(greenland.Lat, greenland.Lon, [1.0 1 1.0])
geoshow([Lat(Selection(1))],[Lon(Selection(1))],'marker','*','color','black','markersize',10)%,'color','black','linewidth',5);

figure(2);
plot(X(Selection)/1000,Y(Selection)/1000, 'k', 'linewidth', 2)
% plot(Lon(Selection),Lat(Selection), 'k', 'linewidth', 2)
title('Traverse Route')
% xlabel('Longitude ^{.\circ}')
% ylabel('Latitude ^{.\circ}')
xlabel('Easting [km]')
ylabel('Northing [km]')
set(gca, 'FontSize', 12, 'FontWeight','bold')
grid on

figure(3);
plot(DistanceKm,Zegm08(Selection),'k','linewidth',2)
xlabel('Distance [km]')
ylabel('WGS84-EGM2008 Elevation [m]','rotation',270)
title('Topographic Profile')
set(gca,'fontsize',14,'fontweight','bold')
grid on

figure(4);
plot3(X(Selection)/1000,Y(Selection)/1000, Zegm08(Selection), 'k','linewidth',2)
title('Topography')
xlabel('Easting [km]')
ylabel('Northing [km]')
zlabel('EGM2008 [mAMSL]')
set(gca, 'FontSize', 12, 'FontWeight','bold')
grid on



