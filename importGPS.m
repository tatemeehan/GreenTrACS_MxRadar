%% importGPS - Loads GreenTrACS GPS Data and Extracs Positions for GPR files.
% This program opens files of three formats. GeoHx and NetR8 are .csv formats
% for 2016 traverse. RTKgps reads .txt format from RTKLIB post-processing.
% Positions are converted from Lat, Lon to PSN. Elevation is corrected from
% WGS84 to Geoid via EGM08. A static shift is applied to correct for antenna 
% height above GrIS surface. Slope, Speed, and Heading are calculated. 
% Data is stored in .mat structure "GeoLocation" to be read into MxRadar.m

%  Boise State University: Tate Meehan GreenTrACS 2017

clear; close all; clc;
addpath './functions'% Load GPS.csv
workingDir = pwd;
% Locate GPS data directory and designate output directory
getDir = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/geolocation';
writeDir = uigetdir(getDir,'Write Directory');
getDir = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/GPS/Processed';
dataDir = uigetdir(getDir,'GPS Data Directory');

isWrite = 1;
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

%% Remove Static GPS Traces

% Select the Removal Method
rmvMethod = questdlg('Which method for static position removal?','Choose Method',...
    'Automatic','Manual','No Removal','Automatic');
% Automatic Removal
if strcmp(rmvMethod,'Automatic')
    win = 11;
    [GPSix,GPSixEdges] = removeStaticGPS(Lon,Lat,win);
end
% Manual Removal
if strcmp(rmvMethod,'Manual')
    
    % Create Figure with Lat and Lon versus Index
    figure(100);clf; hold on;
    plot(1:length(Time),abs(Lat)-mean(abs(Lat)),'k')
    plot(1:length(Time),abs(Lon)-mean(abs(Lon)),'b')
    
    % Allocate
    tmpGPSix = [];
    iter = 1;
    surgIter = 0;
    surgIx = [];
    undo = [];
    
    % Picker Loop is Broken with 'Enter' key
    while 1
        % Manual Pick with Cursor
        [tIx,yIx,b] = ginput(1);
        % Select GPS Edges with Left Click or Space Bar
        if b == 1
            isPick = 1;
        elseif b == 32
            isPick = 1;
        else
            isPick = 0;
        end
        
        % Picking Options:
        % 'Enter' ~ exit the picker
        % 'Left Click' or 'Space Bar' ~ pick edges of good GPS locations
        % 'r' ~ Select edge positions for surgical position removal
        % 'z' ~ 10 seconds of zoom tool
        % 'u' ~ undo last GPS pick
        % 'i' ~ undo last surgical removal pick
        
        if isempty(b)
            break
        elseif ~isPick
            % Neutralize Iteration Count
            iter = iter-1;
            
            % Press 'z' for zoom
            % You have 10 Seconds of Zoom
            if b == 122
                zoom on
                pause(10)
                zoom off
            end
            
            % Surgical Removal of extra stops not associated with file save
            if b == 114
                surgIter = surgIter+1;
                surgIx = [surgIx,round(tIx)];
                figure(100);
                plot(round(tIx),yIx,'db','markersize',7,'linewidth',2)
                undo = [undo,'i'];
            end
            
            % Press 'u' for undo last pick
            if b == 117
                if ~isempty(tmpGPSix)
                % Delete the Pick
                tmpGPSix(iter) = [];
                % Clear the Plotted Point
                h100 = figure(100);
                child = get(gca,'Children');
                undoIx = min(strfind(fliplr(undo),'u'));
                delete(child(undoIx));
                clear('child')                
                % Remove 'u' from undo array
                undoIx = max(strfind(undo,'u'));
                undo(undoIx) = [];
                iter = iter-1;
                end
            end
            
            % Press 'i' for undo last surgical remove index
            if b == 105
                if ~isempty(surgIx)
                % Delete the Pick
                surgIx(surgIter) = [];
                % Clear the Plotted Point
                figure(100);
                child = get(gca,'Children');
                undoIx = min(strfind(fliplr(undo),'i'));                
                delete(child(undoIx));
                clear('child')
                % Remove 'i' from undo array
                undoIx = max(strfind(undo,'i'));
                undo(undoIx) = [];
                surgIter = surgIter-1;
                end
            end
        else
            figure(100); 
            plot(round(tIx),yIx,'rx','markersize',15,'linewidth',2)
            tmpGPSix(iter) = round(tIx);
            undo = [undo,'u'];
        end
        iter = iter+1;
    end
    
    % Check and Sort Surgical Removal 
    if ~isempty(surgIx)
        if mod(length(surgIx),2) ~= 0
            warning('Unpaired Surgery.. will not Operate')
            surgIx = [];
        else
            surgIx = sort(surgIx);
            surgIx = reshape(surgIx,2,length(surgIx)/2);
            nSurg = size(surgIx,2);
            surgIter = 1;
        end

    end
       
    % Check that GPS Segments have Start and Stop Edges
    if mod(length(tmpGPSix),2) == 0
        tmpGPSix = reshape(tmpGPSix,2,length(tmpGPSix)/2);
        GPSix = [];
        GPSixEdges = zeros(size(tmpGPSix));
        for kk = 1:size(tmpGPSix,2)
            % Perform Surgical Position Removal
            if ~isempty(surgIx)
                tmpIx = find(surgIx(1,:)>tmpGPSix(1,kk) & surgIx(2,:)<tmpGPSix(2,kk));
                if isempty(tmpIx)
                    % Append Good GPS Positions
                    GPSix = [GPSix;(tmpGPSix(1,kk):tmpGPSix(2,kk))'];
                    % Catch Dupilcated GPSix
                    GPSix = unique(GPSix);
                else
                    nOps = length(tmpIx);
                    for ll = 1:nOps
                        if (surgIx(1,surgIter) > tmpGPSix(1,kk)) && (surgIx(2,surgIter) < tmpGPSix(2,kk))
                            % Number of Positions to Remove
                            tmp = length(surgIx(1,surgIter)+1:surgIx(2,surgIter)-1);
                            GPSix = [GPSix;([tmpGPSix(1,kk):surgIx(1,surgIter),surgIx(2,surgIter):tmpGPSix(2,kk)])'];
                            % Catch Dupilcated GPSix
                            GPSix = unique(GPSix);
                        end
                        surgIter = surgIter+1;
                    end
                end
                % Find Start/Stop Edges
                GPSixEdges(1,kk) = find(GPSix == tmpGPSix(1,kk));
                GPSixEdges(2,kk) = length(GPSix);
            else
                % Append Good GPS Positions
                GPSix = [GPSix;(tmpGPSix(1,kk):tmpGPSix(2,kk))'];
                GPSix = unique(GPSix);
                % Find Start/Stop Edges
                GPSixEdges(1,kk) = find(GPSix == tmpGPSix(1,kk));
                GPSixEdges(2,kk) = find(GPSix == tmpGPSix(2,kk));
            end
        end
    else
        error('Matched Edges were not Selected. Please Select Again')
    end
    % if manual picker is exited without any picks
    if isempty(tmpGPSix)
        GPSix = 1:length(LLZDT);
        GPSixEdges = [GPSix(1);GPSix(end)];
    end
end
% No Removal
if strcmp(rmvMethod,'No Removal')
    GPSix = 1:length(LLZDT);
    GPSixEdges = [GPSix(1);GPSix(end)];
end

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
R = 5; % Smoothing Window Rank
X = 1:length(Zwgs84); % Estimtes to Smooth

% Smooth Coordinate Estimate
Lon = nonParametricSmooth(X,Lon,X,R);
Lat = nonParametricSmooth(X,Lat,X,R);

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
% % Edit Wagon Wheel Picks : (
% %Temporariy Load the Old Datafile
% % Core7Survey
% % load(['/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/geolocation/5-14-17-Core7Survey','/','5-14-17-Core7Survey-GPS-LINE01_05-tmp.mat'])
% % rmvIx = [[74:82],[354:362],[678:686],[990:1006]];
% % Test
% load(['/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/geolocation/6-3-17-Core12SpurE','/','6-3-17-Core12SpurE-GPS-LINE00_01.mat'])
% load(['/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/geolocation/6-9-17-Core13Core14Traverse','/','6-9-17-Core13Core14Traverse-GPS-LINE01_07.mat'])
% %Core9Spiral
% load(['/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/geolocation/5-17-17-Core9Spiral','/','5-17-17-Core9Spiral-GPS-LINE00_05.mat']);
% Ix = 1:length(GeoLocation.X);
% rmvIx = [[401:411],[564:574],[705:720],[884:895],[1057:1067],[1218:1225]];
% gapIx = Ix;
% gapIx(rmvIx) = [];
% Selection = gapIx;
% tmpGPSixEdges = find(diff(gapIx)>1);
% % Append End
% tmpGPSixEdges = [tmpGPSixEdges,length(gapIx)];
% % Append Top Row
% tmpGPSixEdges = [[1,tmpGPSixEdges(1:end-1)+1];tmpGPSixEdges];


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
%     field12 = 'GPSixEdges';     value12 = tmpGPSixEdges;
%     field13 = 'Slope';          value13 = Slope(Selection);
%     field14 = 'Speed';          value14 = Velocity(Selection);
%     field15 = 'Heading';        value15 = heading(Selection);
%     field16 = 'dX';             value16 = dX(Selection);
%     field17 = 'dY';             value17 = dY(Selection);
%     field18 = 'dZ';             value18 = dZ(Selection);
    field12 = 'GPSixEdges';     value12 = GPSixEdges;
    field13 = 'Slope';          value13 = Slope;
    field14 = 'Speed';          value14 = Velocity;
    field15 = 'Heading';        value15 = heading;
    field16 = 'dX';             value16 = dX;
    field17 = 'dY';             value17 = dY;
    field18 = 'dZ';             value18 = dZ;

% Chunk for WagonWheel Ovewrite
%     field1 = 'X';               value1 = [GeoLocation.X(Selection)];
%     field2 = 'Y';               value2 = [GeoLocation.Y(Selection)];
%     field3 = 'Z_EGM08';         value3 = [GeoLocation.Z_EGM08(Selection)];
% %     field4 = 'Zone';            value4 = {Zone(Selection)};
%     field5 = 'Lon';             value5 = [GeoLocation.Lon(Selection)];
%     field6 = 'Lat';             value6 = [GeoLocation.Lat(Selection)];
%     field7 = 'Z_WGS84';         value7 = [GeoLocation.Z_WGS84(Selection)];
%     field8 = 'N_EGM08';         value8 = [GeoLocation.N_EGM08(Selection)];
%     field9 = 'Date';            value9 = {GeoLocation.Date(Selection)};
%     field10 = 'Time';           value10 = {GeoLocation.Time(Selection)};
% %     field11 = 'GPSix';          value11 = GPSix(Selection);
%     field12 = 'GPSixEdges';     value12 = tmpGPSixEdges;
%     field13 = 'Slope';          value13 = [GeoLocation.Slope(Selection)];
%     field14 = 'Speed';          value14 = [GeoLocation.Speed(Selection)];
%     field15 = 'Heading';        value15 = [GeoLocation.Heading(Selection)];
%     field16 = 'dX';             value16 = [GeoLocation.dX(Selection)];
%     field17 = 'dY';             value17 = [GeoLocation.dY(Selection)];
%     field18 = 'dZ';             value18 = [GeoLocation.dZ(Selection)];
    
    GeoLocation = struct(field1,value1,field2,value2,field3,value3,...%field4,value4,...
        field16,value16,field17,value17,field18,value18,...
        field5,value5,field6,value6,field7,value7,field8,value8,...
        field13,value13,field14,value14,field15,value15,...
        field9,value9,field10,value10,...
        field12,value12);%,field11,value11);

% Save GPS Data
if isWrite
cd(writeDir)
putName = input('Input File Name (Exclude the File Extension): ','s');
writeName = [putName,'.mat'];
% writeName = 'MxRadarGPSCore15SpurW061317.mat';
save(writeName,'GeoLocation','-v7.3')
cd(workingDir)
end

%% Create Maps

figure(1);clf;
ax = worldmap('Greenland');
greenland = shaperead('landareas', 'UseGeoCoords', true,...
  'Selector',{@(name) strcmp(name,'Greenland'), 'Name'});
patchm(greenland.Lat, greenland.Lon, [1.0 1 1.0])
geoshow([Lat(Selection(1))],[Lon(Selection(1))],'marker','*','color','black','markersize',10)%,'color','black','linewidth',5);

figure(2);clf
% plot(X(Selection)/1000,Y(Selection)/1000, 'k', 'linewidth', 2)
% plot(X/1000,Y/1000, 'k', 'linewidth', 2)
plot(GeoLocation.X/1000,GeoLocation.Y/1000, '.k', 'linewidth', 2)

% plot(Lon(Selection),Lat(Selection), 'k', 'linewidth', 2)
title('Traverse Route')
% xlabel('Longitude ^{.\circ}')
% ylabel('Latitude ^{.\circ}')
xlabel('Easting [km]')
ylabel('Northing [km]')
set(gca, 'FontSize', 12, 'FontWeight','bold')
grid on
hold on;
for kk = 1:size(GeoLocation.GPSixEdges,2)
%     plot(X(GPSixEdges(1,kk):GPSixEdges(2,kk))/1000,Y(GPSixEdges(1,kk):GPSixEdges(2,kk))/1000, 'linewidth', 2)
    plot(GeoLocation.X(GeoLocation.GPSixEdges(1,kk):GeoLocation.GPSixEdges(2,kk))/1000,GeoLocation.Y(GeoLocation.GPSixEdges(1,kk):GeoLocation.GPSixEdges(2,kk))/1000, 'linewidth', 2)

end


figure(3);clf
plot(DistanceKm,Zegm08(Selection),'k','linewidth',2)
xlabel('Distance [km]')
ylabel('WGS84-EGM2008 Elevation [m]','rotation',270)
title('Topographic Profile')
set(gca,'fontsize',14,'fontweight','bold')
grid on

figure(4);clf
plot3(X(Selection)/1000,Y(Selection)/1000, Zegm08(Selection), 'k','linewidth',2)
title('Topography')
xlabel('Easting [km]')
ylabel('Northing [km]')
zlabel('EGM2008 [mAMSL]')
set(gca, 'FontSize', 12, 'FontWeight','bold')
grid on



