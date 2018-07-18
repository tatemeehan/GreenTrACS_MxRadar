clear; %close all; clc;
%% MultiOffsetRoutine Reads, Sorts, and Processes the Multiplexed Sensors
%  and Software MultiChannel GPR Record. Processed data are stored in the
%  common-offset domain. Offsets are sorted into common shot gathers for
%  Surface Velocity Analysis performed by Monte Carlo BootStrapping and 
%  Iteratively ReWeighted Least-Squares Method in offset-traveltime space. 
%   
%  This Code is Parallelized - User Must Specify Number of Workers
%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%
%  Written by Tate Meehan, Boise State University, GreenTrACS 2016-2017   %  
%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%

% Input  Data Directory

% 500 MHz Data
% GreenTrACS 2016
% dataDir = '/home/tatemeehan/GreenTracs2016/GPR_DATA/PulseEKKO/500MHz/6-4-16-Core7-Spiral';
dataDir = '/home/tatemeehan/GreenTracs2016/GPR_DATA/PulseEKKO/500MHz/6-2-16-Core7-Spur-W';

% GreenTrACS 2017
% dataDir = '/home/tatemeehan/GreenTracs2017/6-19-17-SummitRouteII';
% dataDir = 'E:\GreenTrACS\2017\6-19-17-SummitRouteII';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\6-19-17-SummitRouteII';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\6-18-17-SummitRouteI';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\6-15-17-Core15Core16Traverse';
% dataDir = '/home/tatemeehan/GreenTracs2017/6-13-17-Core15SpurW';
% dataDir = '/sonichome/tatemeehan/GreenTracs2017/6-13-17-Core15SpurW';
% dataDir = 'E:\GreenTrACS\2017\6-13-17-Core15SpurW';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\6-13-17-Core15SpurW';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\6-12-17-Core14Core15Traverse';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\6-10-17-Core14SpurW';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\6-9-17-Core13Core14Traverse';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\6-8-17-Core13SpurE';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\6-7-17-Core13SpurW';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\6-6-17-Core12Core13Traverse';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\6-3-17-Core12SpurE';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\6-2-17-Core12SpurW';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\6-1-17-Core11Core12Traverse';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\5-31-17-Core11SpurW';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\5-27-17-Core11SpurE';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\5-26-17-Core11Spiral';
% dataDir = '/run/media/tatemeehan/ONE/GreenTrACS2017/PulseEKKO/500MHz/5-26-17-Core11Spiral';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\5-24-17-Core10SpurE';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\5-25-17-Core10SpurW';
% dataDir = '/run/media/tatemeehan/ONE/GreenTrACS2017/PulseEKKO/500MHz/5-26-17-Core10Core11Traverse';
% dataDir = '/home/tatemeehan/GreenTracs2017/5-26-17-Core10Core11Traverse';
% dataDir = 'D:\GreenTrACS\2017\GPR_DATA\PulseEKKO\500MHz\5-17-17-Core9Spiral';
% dataDir = 'D:\GreenTrACS\2017\GPR_DATA\PulseEKKO\500MHz\5-11-17-Core8Test';
% dataDir = 'D:\GreenTrACS\2017\GPR_DATA\PulseEKKO\500MHz\5-4-17-SummitRun';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\5-5-17-SummitCal';
% dataDir = '/home/tatemeehan/GreenTracs2017/GPR_DATA/PulseEKKO/500MHz/5-5-17-SummitCal';
% dataDir = 'E:\GreenTrACS\GPR_DATA/PulseEKKO/500MHz/6-4-16-Core7-Spiral';
% dataDir = '/home/tatemeehan/GreenTracs2017/GPR_DATA/PulseEKKO/500MHz/5-17-17-Core9Spiral';

% 1GHz Data
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\1GHz\6-19-17-SummiteRouteII';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\1GHz\6-7-17-Core13SpurW';
% dataDir = 'D:\GreenTrACS\2017\GPR_DATA\PulseEKKO\1GHz\5-17-17-Core9Spiral';
% dataDir = 'D:\GreenTrACS\2017\GPR_DATA\PulseEKKO\1GHz\5-11-17-Core8Test';
% dataDir = 'D:\GreenTrACS\2017\GPR_DATA\PulseEKKO\1GHz\5-4-17-SummitRun';
% dataDir = 'D:\GreenTrACS\2017\GPR_DATA\PulseEKKO\1GHz\5-5-17-SummitCal';
% dataDir = '/home/tatemeehan/GreenTracs2017/GPR_DATA/PulseEKKO/1GHz/5-10-17-Core8SpurE';

% Add additional useful pathways
% addpath 'E:\GreenTrACS\2017\GPR_Processing\MultiOffset\Save';
% addpath 'E:\GreenTrACS\2017\GPR_Processing\MultiOffset\SeismicLab\SeismicLab\codes\fx';
% addpath 'D:\GreenTrACS\2017\GPR_Processing\MultiOffset\Save';
% addpath '/run/media/tatemeehan/RED/GreenTrACS/2017/GPR_Processing/MultiOffset/Save'
addpath '/home/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/Save'
addpath '/home/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/TM'
% addpath '/home/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/SeismicLab/SeismicLab/codes/fx'
% addpath '/sonichome/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/Save'
% addpath '/sonichome/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/SeismicLab/SeismicLab/codes/fx'

%% Control Floe
% Parallel Computing Enabled
isParallel = 1;

% Read Data
isReadSensorsSoftware = 1;     % Read Multiplexed Data
isLoadTimeHorizons = 0;        % Load Previously Picked Time Horizons

% Export Data
isWriteTimeHorizons = 0;

% Process Data
isTrimTWT = 0;          % Truncate Recorded Data for Near-Surface Analysis
isKill = 0;             % Kill Unanted Channels
isFXdecon = 1;          % Fx-Predicitive Deconvolution
isSWEDISH = 1;          % Perform Surface Velocity Analysis

%% Control Parallelization
if isParallel
    % Wake Parallel Computing
    if isempty(gcp('nocreate'))     
        nWorkers = 9;
        p = parpool(nWorkers);
    else nWorkers = 9;
    end
else
    nWorkers = 1;
end

%% Import Meta Data
load('CalibrationChannelShiftTrough500MHz2017.mat');
chanShift = Calibration(1).chanShift; % Import chanShift
load('LateNite.mat');   % Import Cool Colormap

lineNo = [0,1];           % Array of data "LINE" numbers
nFiles = length(lineNo);    % Number of Files
nChan = 9;             % Number of Recorded Channels
chan =  1:nChan;       % Linear Array of Record Channels
liveChan = chan;

% Establish Tx Rx Geometry for CMP Gathering
nTx = 3;               % Number of Transmitters in Sequence
nRx = 3;               % Number of Receivers in Sequence

% Allocate Memory
time = cell(1,nFiles);
%% Import GPS Information
% Yet To be Completed
%% Read GPR Data
if isReadSensorsSoftware
    Rad = cell(1,nFiles);
    trhd = cell(1,nFiles);
    for ii = 1 : nFiles
        tic
        %------------------------------------------------------------------
        % Multiplexed Channel Record
        if lineNo(ii)<10
            filename=['LINE0' num2str(lineNo(ii))];
        else
            filename=['LINE' num2str(lineNo(ii))];
        end
        
        % Read Data
        filepath = fullfile(dataDir,filename);
        [Rad{ii},hdr1,trhd{ii},dt,f0,~,dx] = readSensorsSoftwareData( filepath );
        
        [~, multiplexNtrcs] = size(Rad{ii});% No. Traces in Multiplexed Data
        
        % Install Transmitter and Receiver Sequencing and Geometry
        % 500 MHz Offsets
        if f0 == 500
            if nChan == 9
                % 9 Channel Offsets
                txGeo = [0, 0, 0, 1.33, 1.33, 1.33, 2.67, 2.67, 2.67]; % Tx Sequence & Absolute Position
                rxGeo = [4, 8, 12, 4, 8, 12, 4, 8, 12]; % Rx Sequence & Absolute Position
            elseif nChan == 6
                % 6 Channel Offsets
                txGeo = [1.33, 1.33, 1.33, 2.67, 2.67, 2.67]; % Tx Sequence & Absolute Position
                rxGeo = [4, 8, 12, 4, 8, 12]; % Rx Sequence & Absolute Position
            end
            % 1 GHz Offsets
        elseif f0 == 1000
            txGeo = [0, 0, 0, 0.67, 0.67, 0.67, 1.33, 1.33, 1.33]; % Tx Sequence & Absolute Position
            rxGeo = [2, 4, 6, 2, 4, 6, 2, 4, 6]; % Rx Sequence & Absolute Position
        end
            
        offsetArray = rxGeo - txGeo;            % Offset Array for CMP Gather
        nChan = length(offsetArray);            % Refresh nChan before kill
        farOffset = max(offsetArray);           % Determine Far Offset
        farChan = find(offsetArray == farOffset);
        nearOffset = min(offsetArray);          % Determine Near Offset
        nearChan = find(offsetArray == nearOffset);
        
        % Determine Data Acquisition  Method
        if all(trhd{ii}(5,:) == 0)
            isFreeRun = 1;
            isOdometer = 0;
        else
            isOdometer = 1;
            isFreeRun = 0;        
        end
        
        % Remove Static Traces if is Free Run Acquisition
        if isFreeRun
            % Process Near Offset Channel
            disp(' ')
            fprintf('Begin Static Trace Removal \n')
            tic
            nearRad = processTraceRemoval... % Processes Near Offset Data
                (Rad{ii}(:,nearChan:nChan:end), f0, dt );
            % Remove Duplicate Static Traces
            dupIx = removeStaticTrace( nearRad, multiplexNtrcs, nearChan, nChan );
            
            trhd{ii}(:,dupIx) = [];         % Remove Static Trace Headers from Multiplexed Record
            trhd{ii}(1,:) = 1:length(trhd{ii}); % Configure Trace Indicies
            trhd{ii}(2,:) = [0:dx:(length(trhd{ii})-1)*dx]; % Configure Distance
            Rad{ii}(:,dupIx) = [];      % Remove Static Traces from Multiplexed Data
            xArray = trhd{ii}(2,:);         % Define Configured Distance xArray
            
            fprintf('Static Trace Removal Done \n')
            toc
            display(' ')
            
        end
        
        % Remove Skip Traces if Wheel Odometer Acquisition
        if isOdometer
        dupIx = find(~(diff(trhd{ii}(23,:)))); % Find Skipped Traces
        
        for jj = dupIx
            trhd{ii}(1,jj:end) = trhd{ii}(1,jj:end) - 1; % ReConfigure Trace Indicies
%             trhd{ii}(2,jj:end) = trhd{ii}(2,jj:end) - dx;% ReConfigure Distance
        end
        
        trhd{ii}(:,dupIx) = []; % Remove Skipped Traces from Trace Header
        Rad{ii}(:,dupIx) = []; % Remove Skipped Traces from Multiplexed Data
        xArray = trhd{ii}(2,:); % Define ReConfigured Distance as xArray
        end
        % Gain Far Channels for Cable Attenuation
        if f0 == 500
            gainIx = find(trhd{ii}(23,:) == 3 | trhd{ii}(23,:) == 6 | trhd{ii}(23,:) == 9);
            Rad{ii}(:,gainIx) = Rad{ii}(:,gainIx).*1.7;
        end
                
        % Kill Unwanted Channels
        if isKill
            if ii == 1;
                killChan = [7];         % Kill Channel Number
                liveChan(killChan) = [];% Array of Live Channels
                % Update Killed Channel Shifts
                chanShift(killChan) = [];
%                 nWorkers = length(chan);    %re-assign parallel Workers
            end
            killIx = find(ismember(trhd{ii}(23,:),killChan));
            % Kill Traces
            Rad{ii}(:,killIx) = [];  
            % Remove Killed Trace Headers
            trhd{ii}(:,killIx) = [];        
            % Configure Trace Indicies
            trhd{ii}(1,:) = 1:length(trhd{ii}); 
             % Define Configured Distance xArray
            xArray = trhd{ii}(2,:);        
            % Remove Killed Offsets
            offsetArray(killChan) = []; 
            % Number of Live Channels
            nChan = nChan - length(killChan);
            % Looping Array of Record Channels
            chan =  1:nChan;
        end
        
        % Pad Data with Instrument Zero
        padding = 100;
        instrumentPad = zeros(padding,size(Rad{ii},2));
        if padding ~= 0
            for jj = 1:size(Rad{ii},2)
                instrumentZero = Rad{ii}(1,jj);
                instrumentPad(:,jj) = ones(padding,1).*instrumentZero;
            end
            Rad{ii} = [instrumentPad;Rad{ii}];
        end  
        
        % Trim Time Window
        if isTrimTWT
            reSample = 450 + padding;   % Number of Wanted Samples
            Rad{ii} = Rad{ii}(1:reSample,:);
        end
        
        % Allocation Here
        if ii == 1;
            Radar = cell(nChan,nFiles); traceIx = cell(nChan,nFiles);
            Array = cell(nChan,nFiles);
        end
        
        parfor (jj =  1:nChan, nWorkers)
%         for jj = chan
            % DeMux Sequential Data
            [Radar{jj,ii}, traceIx{jj,ii}, Array{jj,ii}] = DeMux(Rad{ii}, trhd{ii}, liveChan(jj));
            
            % Extract Full-fold Traces & Sort Antenna Positions
            gatherLength = length(Array{jj,ii}(1,:)); % Length of Each Chan
            traceMod = mod(length(xArray),nChan); % Un-folded Channel Index
            
            % Flag Un-binned Traces
            if jj <= traceMod
                xTrc = 1;
            else
                xTrc = 0;
            end
                
                % Remove un-Binned Common Offset Traces
                Radar{jj,ii} = Radar{jj,ii}(:,1:gatherLength - xTrc);
                
                % Create Temporary Two-Way Time
                tmpTime = [0:dt:(size(Radar{jj,ii},1)-(1+padding))*dt];

                % Apply Systematic Channel Shifts
                if chanShift(jj) < 0
                    tmp = padarray(Radar{jj,ii},[abs(chanShift(jj)),0],...
                        mean(Radar{jj,ii}(1,:)),'pre');
                    Radar{jj, ii} = tmp(1:length(tmpTime),:);

                elseif chanShift(jj) > 0
                    tmp = padarray(Radar{jj,ii},[abs(chanShift(jj)),0],...
                        mean(Radar{jj,ii}(1,:)),'post');
                    Radar{jj,ii} = tmp(abs(chanShift(jj))+1:end-padding,:);
                   
                else
                    Radar{jj,ii} = Radar{jj,ii};
                  
                end
                
                % Process Common Offset Channels
                disp(' ')           
                fprintf(['Begin Signal Processing in Common-Offset Domain ',...
                    filename, ' CHAN0', num2str(jj),'\n'])

                Radar{jj,ii} = processCommonOffset(Radar{jj,ii}, f0, dt );

        end
        fprintf('Signal Processing Done \n')
        toc
        display(' ')
    end
end
%% Concatenate Files
    % Grab Early Time Data For AirWave Cross Correlation Alignment
    % 2016
        xcorrWindow = [30, 105, 180, 0, 80, 145, -20, 45, 130;...
                    60, 150, 210, 40, 110, 175, 5, 80, 160] + padding;
    % 2017
%     xcorrWindow = [175, 250, 85, 140, 230, 75, 100, 180, 25;...
%                     200, 300, 150, 160, 260, 115, 150, 220, 75] + padding;
    if isKill
        xcorrWindow(:,killChan) = [];
    end
    
    % Cross-Correlation for multipleFile Alignment Prior to Concatenation
    Radar = xcorrAlignTraces( Radar, xcorrWindow );
    
    % Concatenate Files
    isCat = 1;
    if isCat
        nFiles = 1;
        catRadar = cell(nChan,nFiles);
        for jj = chan            
            catRadar{jj} = [Radar{jj,:}];
        end
         Radar = catRadar; clear catRadar;
    end
    
%% Branch Processing for Direct and Reflected Arrivals
directRadar = Radar;

%% Median Background Subtraction Filter for Reflection Analysis
    medR = 1000; % Rank of Median Subtraction Window: W = 2R+1
    for ii = 1:nFiles
        parfor (jj = chan, nWorkers)
%         for jj = chan;
            Radar{jj,ii} = movingMedianSubtraction(Radar{jj,ii},medR);
%             pickRadar{jj,ii} = pickRadar{jj,ii}-median(pickRadar{jj,ii},2)...
%                 *ones(1,size(pickRadar{jj,ii},2));
        end
    end
    
%% FX Predicive Deconvolution
if isFXdecon
    fprintf('Begin FX-Predictive Deconvolution \n')
    tic
    
    % Select Temporal or Spatial Deconvolution Gates: 1 on; 0 off;
    isTemporalFXdecon = 1;
    isSpatialFXdecon = 0;

    FXdeconTX
    
    fprintf('FX-Deconvolution Done \n')
    toc
    display(' ')
end

%% Snow Water Equivalent and Density for Ice Sheet Height
if isSWEDISH
    if isLoadTimeHorizons
        TimeHorizonFilename = '';
        load(TimeHorizonFilename);
        % Determine Number of Direct Wave Horizons
        nDirectHorizon = size(DirectFBpick,2);
        % Determine Number of Reflection Horizons
        nReflectionHorizon = size(ReflectionFBpick,2);
    else
    %% Semi-Automatic Radar Wave Picking
    fprintf('Begin PolarPicker \n')
    display(' ')
    tic

    % AGC Gain for PolarPicker 
    for ii = 1:nFiles
        parfor (jj = chan, nWorkers)
            directRadar{jj,ii} = AGCgain(directRadar{jj,ii},50,2);
        end
    end
    
    % Pick Direct Wave Arrival
    display('Pick Direct Wave')
    display(' ')
    display('Initial Horizon must be the Air Wave!')
    display(' ')
    [~, DirectFBpick] = polarPickerT8(directRadar);
    
    reflectRadar = Radar;
    
        % AGC Gain for PolarPicker 
    for ii = 1:nFiles
        parfor (jj = chan, nWorkers)
            reflectRadar{jj,ii} = AGCgain(reflectRadar{jj,ii},250,2);
        end
    end
    
    % Plot Common-offset gathers
    isPlotOffsetGathers = 1;
    if isPlotOffsetGathers
        for ii = 1:nFiles
            for jj = chan
                figure(chan(jj)+(ii-1).*length(chan));...
                    imagesc(reflectRadar{jj,ii});colormap(LateNite);
            end
        end
        
        % Pause to Examine Common-offset Gathers
        display('Review the Common-Offet Gathers')
        display('Press F5 to Continue')
        display(' ')
        keyboard
    end
    
    % Pick Primary Reflection Arrival
    display('Pick Primary Reflection')
    display(' ')
    display('Reflection Horizons Must be Picked Sequenctially in Time!')
    display(' ')
    [~, ReflectionFBpick] = polarPickerT8(reflectRadar);
    
    clear reflectRadar
    
    % Determine Number of Direct Wave Horizons
    nDirectHorizon = size(DirectFBpick,2);
    % Determine Number of Reflection Horizons
    nReflectionHorizon = size(ReflectionFBpick,2);
        
    % Convert Samples to ns
    for ii = 1:nFiles
        time{ii} = [0:dt:(size(Radar{ii,1},1)-(1+padding(ii)))*dt];
        for jj = chan
            for dh = 1:nDirectHorizon
                DirectFBpick{jj,dh,ii} = DirectFBpick{jj,dh,ii}.*dt;
            end
            for rh = 1:nReflectionHorizon
                ReflectionFBpick{jj,rh,ii} = ReflectionFBpick{jj,rh,ii}.*dt;
            end
        end
    end
    
    if isWriteTimeHorizons
        save('6-2-16-Core7-Spur-W-TimeHorizon.mat','DirectFBpick','ReflectionFBpick','-v7.3');
    end

    fprintf('PolarPicker Done \n')
    toc
    display(' ')
    end
%% Perform Horizon Velocity Analysis for Estimates of Density, Depth, & SWE
    fprintf('Begin Horizon Velocity Analysis \n')
    display(' ')
    tic
% Run Memory Allocation    
HVAmemoryAllocation

% Extract AirWave Picks and Perform Residual Subtraction Velocity Analysis
%     parfor (ii = 1:nFiles, nWorkers - (nWorkers-nFiles)) 
    for ii = 1:nFiles
        looper = 1:size(Radar{ii},2);
        % Concatenate Air Wave Picks
        for dh = 1:nDirectHorizon
            % Concatenate Primary Reflection Picks for Horizon hh
            GatherDirectPicks{dh,ii} = cat(2,DirectFBpick{:,dh,ii});
            DirectPicks = GatherDirectPicks{dh,ii};
            % Air Wave Velocity Analysis
            if dh  == 1;
%                 for jj = looper
                 parfor (jj = looper, nWorkers)
%                             for jj = 1:length(Radar{ii})
                    
                    % Extract Picks
                    AirPick = DirectPicks(jj,:);
                    
                    for kk = 1:250 % 250 Random Draws
                        % Cross-Validation for Surface Velocity Estimation 7-28-17
                        nCut = randsample([0,1,2],1);
                        cutChan = randsample(liveChan,nCut);
                        xvalChan = liveChan;
                        cutIx = find(ismember(liveChan,cutChan));
                        xvalChan(cutIx) = [];
                        xvalIx = find(ismember(liveChan,xvalChan));
                        
                        xvalAirPick = AirPick(xvalIx);
                        xvalOffset = offsetArray(xvalIx);
                        
                        % Compute Air Wave Arrival Velocity and Intercept Time
                        % OLS Scheme
                        % [xVair{ii,jj}(kk,1), ~] = DirectWave(xvalOffset,xvalAirPick);
                        % IRLS Scheme
                        [xVdir{dh,jj,ii}(kk,1), xToDir{dh,jj,ii}(kk,1)] = DirectWaveIrls(xvalOffset,xvalAirPick);
                    end
                    
                    % Residual Error Analysis
                    VoAir = mean([xVdir{dh,jj,ii}]);
                    
                    % Find Velocity Residual for Additional Trace Shifts & Calulate Residual Time
                    deltaT{ii,jj} = (offsetArray./VoAir) - (offsetArray./0.299);
                    ResAirPick = AirPick - deltaT{ii,jj};
                    % Recalculate Airwave Moveout Velocity
                    [VoDir{dh,jj,ii}, toDir{dh,jj,ii}] = DirectWaveIrls(offsetArray,ResAirPick);

                    %                 xAirTo{ii,jj} = xToDir{dh,jj,ii};
                    AirTo{ii,jj} = toDir{dh,jj,ii};
                    
                    
                    % Apply Residual Time Shifts
                    %                 ResGatherDirPicks{dh,jj,ii} = DirectPicks(jj,:) - deltaT{ii,jj};
                    
                    % Estimate Wavelet Depth
                    %             xDepthDir{dh,jj,ii} = xVdir{dh,jj,ii}.*(xToDir{dh,jj,ii}-xAirTo{ii,jj});
                    xDepthDir{dh,jj,ii} = 0;
                    % Estimate Direct Wave Density
                    %             xRhoDir{dh,jj,ii} = DryCrim([xVdir{dh,jj,ii}]);
                    xRhoDir{dh,jj,ii} = 0;
                    % Error Analysis
                    %             VoDir{dh,jj,ii} = mean([xVdir{dh,jj,ii}]);
                    VoVarDir{dh,jj,ii} = var([xVdir{dh,jj,ii}]);
                    %             toDir{dh,jj,ii} = mean([xToDir{dh,jj,ii}]);
                    toVarDir{dh,jj,ii} = var([xToDir{dh,jj,ii}]);
                    %             DepthDir{dh,jj,ii} = mean([xDepthDir{dh,jj,ii}]);
                    DepthDir{dh,jj,ii} = VoDir{dh,jj,ii}.*toDir{dh,jj,ii};
                    %             RhoDir{dh,jj,ii} = mean([xRhoDir{dh,jj,ii}]);
                    RhoDir{dh,jj,ii} = DryCrim(VoDir{dh,jj,ii});
                    %             CovZP = cov([xDepthDir{dh,jj,ii}],[xRhoDir{dh,jj,ii}]);
                    CovZP = zeros(2,2);
                    DepthVarDir{dh,jj,ii} = CovZP(1,1);
                    RhoDirVar{dh,jj,ii} = CovZP(2,2);
                    CovDepthRhoDir{dh,jj,ii} = CovZP(1,2);
                end
            else
                % Cross-Validation for Surface Velocity Estimation
%                 for jj = looper
                 parfor (jj = looper, nWorkers)
                    % Multiple Shot Gathers in Population
                    isManyShots = 1;
                    if isManyShots
                        shotRange = 5;
                        ranger = sqrt(([1:size(Radar{ii},2)]-jj).^2); % Compute Distance
                        getIx = find(ranger<=shotRange); % Find Nearby Picks
                        DirPick = DirectPicks(getIx,:) - vertcat(deltaT{ii,getIx}); % Residual Static Corection
                        DirPick = DirPick - vertcat(AirTo{ii,getIx})*ones(1,nChan); % Time-Zero Static Shift
                        pickPool = size(DirPick,1);
                        xvalPool = 1:pickPool;
                        % Cross-Validation for Surface Velocity Estimation 10-12-17
                        for kk = 1:1250 % 1250 Random Draws
                            nCut = randsample([0,1,2],1);
                            cutChan = randsample(liveChan,nCut);
                            xvalChan = liveChan;
                            cutIx = find(ismember(liveChan,cutChan));
                            xvalChan(cutIx) = [];
                            xvalIx = find(ismember(liveChan,xvalChan));
                            pickPool = length(xvalIx);
                            xvalBin = randsample(xvalPool,pickPool,'true');
                            xvalDirPick = ...
                                DirPick(sub2ind(size(DirPick),xvalBin,xvalIx));
                            xvalOffset = offsetArray(xvalIx);
                            
                            % Compute Air Wave Arrival Velocity and Intercept Time
                            % OLS Scheme
                            % [xVair{ii,jj}(kk,1), xToAir{ii,jj}(kk,1)]...
                            %   = DirectWave(xvalOffset,xvalAirPick);
                            % IRLS Scheme
                            [xVdir{dh,jj,ii}(kk,1), xToDir{dh,jj,ii}(kk,1)]...
                                = DirectWaveIrls(xvalOffset,xvalDirPick);
                        end
                    end
                    % Single Shot Gather in Population
                    isSingleShot = 0;
                    if isSingleShot
                        % Re-Extract Direct Wave Picks after Residual Shifts
                        DirPick = DirectPicks(jj,:) - deltaT{ii,jj}; % Residual Static Corection
                        DirPick = DirPick - AirTo{ii,jj}*ones(1,nChan); % Time-Zero Static Shift                       
                        for kk = 1:250 % 250 Random Draws
                            % Cross-Validation for Surface Velocity Estimation 7-28-17
                            nCut = randsample([0,1,2],1);
                            cutChan = randsample(liveChan,nCut);
                            xvalChan = liveChan;
                            cutIx = find(ismember(liveChan,cutChan));
                            xvalChan(cutIx) = [];
                            xvalIx = find(ismember(liveChan,xvalChan));
                            
                            xvalDirPick = DirPick(xvalIx);
                            xvalOffset = offsetArray(xvalIx);
                            
                            % Compute Air Wave Arrival Velocity and Intercept Time
                            % OLS Scheme
                            % [xVair{ii,jj}(kk,1), xToAir{ii,jj}(kk,1)]...
                            %   = DirectWave(xvalOffset,xvalAirPick);
                            % IRLS Scheme
                            [xVdir{dh,jj,ii}(kk,1), xToDir{dh,jj,ii}(kk,1)]...
                                = DirectWaveIrls(xvalOffset,xvalDirPick);
                        end
                    end
                    
                    % Estimate Wavelet Depth
                    xDepthDir{dh,jj,ii} = xVdir{dh,jj,ii}.*abs(xToDir{dh,jj,ii});
                    % Estimate Direct Wave Density
                    xRhoDir{dh,jj,ii} = DryCrim([xVdir{dh,jj,ii}]);
                    % Error Analysis
                    VoDir{dh,jj,ii} = mean([xVdir{dh,jj,ii}]);
                    VoVarDir{dh,jj,ii} = var([xVdir{dh,jj,ii}]);
                    toDir{dh,jj,ii} = mean(abs([xToDir{dh,jj,ii}]));
                    toVarDir{dh,jj,ii} = var([xToDir{dh,jj,ii}]);
                    DepthDir{dh,jj,ii} = mean([xDepthDir{dh,jj,ii}]);
                    RhoDir{dh,jj,ii} = mean([xRhoDir{dh,jj,ii}]);
                    CovZP = cov([xDepthDir{dh,jj,ii}],[xRhoDir{dh,jj,ii}]);
                    DepthVarDir{dh,jj,ii} = CovZP(1,1);
                    RhoDirVar{dh,jj,ii} = CovZP(2,2);
                    CovDepthRhoDir{dh,jj,ii} = CovZP(1,2);
                end
                
            end
            
            % Concatenate Direct Travel Time
            DirectTo{dh,ii} = [toDir{dh,:,ii}];
            DirectToVar{dh,ii} = [toVarDir{dh,:,ii}];
            % Concatenate Direct Velocity
            DirectVelocity{dh,ii} = [VoDir{dh,:,ii}];
            DirectVelocityVar{dh,ii} = [VoVarDir{dh,:,ii}];
            % Concatenate Direct Wave Depth
            DirectDepth{dh,ii} = [DepthDir{dh,:,ii}];
            DirectDepthVar{dh,ii} = [DepthVarDir{dh,:,ii}];
            % Concatenate Reflection Density
            DirectDensity{dh,ii} = [RhoDir{dh,:,ii}];
            DirectDensityVar{dh,ii} = [RhoDirVar{dh,:,ii}];
            % Concatenate Covariance
            CovarianceDepthDensityDirect{dh,ii} = [CovDepthRhoDir{dh,:,ii}];
            
            % Estimte DirectSWE
            DirectSWE{dh,ii} = DirectDepth{dh,ii}.*DirectDensity{dh,ii};
            % Error Propagation Equation
            DirectSWEvar{dh,ii} = DirectDepthVar{dh,ii}.*DirectDensity{dh,ii}.^2 ...
                + DirectDensityVar{dh,ii}.*DirectDepth{dh,ii}.^2 ...
                + 2.*DirectDepth{dh,ii}.*DirectDensity{dh,ii}.*CovarianceDepthDensityDirect{dh,ii};
            
            % Smooth Esitmates
            smoothR = 251;
            dhTWT{dh,ii} = nonParametricSmooth( 1:length(DirectTo{dh,ii}),...
                DirectTo{dh,ii},1:length(DirectTo{dh,ii}),smoothR);
            dhTWTvar{dh,ii} = nonParametricSmooth( 1:length(DirectToVar{dh,ii}),...
                DirectToVar{dh,ii},1:length(DirectToVar{dh,ii}),smoothR);
            dhSnowWaterEqv{dh,ii} = nonParametricSmooth( 1:length(DirectSWE{dh,ii}),DirectSWE{dh,ii},...
                1:length(DirectSWE{dh,ii}),smoothR);
            dhSnowWaterEqvVar{dh,ii} = nonParametricSmooth( 1:length(DirectSWEvar{dh,ii}),...
                DirectSWEvar{dh,ii},1:length(DirectSWEvar{dh,ii}),smoothR);
            dhDensity{dh,ii} = nonParametricSmooth( 1:length(DirectDensity{dh,ii}),...
                DirectDensity{dh,ii},1:length(DirectDensity{dh,ii}),smoothR);
            dhDensityVar{dh,ii} = nonParametricSmooth( 1:length(DirectDensityVar{dh,ii}),...
                DirectDensityVar{dh,ii},1:length(DirectDensityVar{dh,ii}),smoothR);
            dhDepth{dh,ii} = nonParametricSmooth( 1:length(DirectDepth{dh,ii}),...
                DirectDepth{dh,ii},1:length(DirectDepth{dh,ii}),smoothR);
            dhDepthVar{dh,ii} = nonParametricSmooth( 1:length(DirectDepthVar{dh,ii}),...
                DirectDepthVar{dh,ii},1:length(DirectDepthVar{dh,ii}),smoothR);
        end
    end

% Primary Reflection Velocity Analysis               
    for ii = 1:nFiles
        looper = 1:size(Radar{ii},2);
        for rh = 1:nReflectionHorizon
        % Concatenate Primary Reflection Picks for Horizon hh
        GatherReflectionPicks{rh,ii} = cat(2,ReflectionFBpick{:,rh,ii});
        Reflection = GatherReflectionPicks{rh,ii};
%         for jj = looper
        parfor (jj = looper, nWorkers)
            % Jackknife Simulation for Reflection Velocity Estimation
            % Multiple Shot Gathers in Population
            isManyShots = 1;
            if isManyShots
                shotRange = 5;
                ranger = sqrt(([1:size(Radar{ii},2)]-jj).^2); % Compute Distance
                getIx = find(ranger<=shotRange); % Find Nearby Picks
                RefPick = Reflection(getIx,:) - vertcat(deltaT{ii,getIx}); % Residual Static Corection
                RefPick = RefPick - vertcat(AirTo{ii,getIx})*ones(1,nChan); % Time-Zero Static Shift
                pickPool = size(RefPick,1);
                xvalPool = 1:pickPool;
%                 stdRefPick = std(RefPick,0,1);
            end
            % Single Shot Gather in Population
            isSingleShot = 0;
            if isSingleShot
                RefPick = Reflection(jj,:) - deltaT{ii,jj}; % Residual
                RefPick = RefPick - AirTo{ii,jj};           % Bulk Shift
            end
                % Jackknife Simulation for Reflection Velocity Estimation 10-12-17
                if isManyShots
                    for kk = 1:1250 % 1250 Random Draws
                    nCut = randsample([0,1,2],1);
                    cutChan = randsample(liveChan,nCut);
                    xvalChan = liveChan;
                    cutIx = find(ismember(liveChan,cutChan));
                    xvalChan(cutIx) = [];
                    xvalIx = find(ismember(liveChan,xvalChan));
                    pickPool = length(xvalIx);
                    xvalBin = randsample(xvalPool,pickPool,'true');
                    xvalRefPick = ...
                        RefPick(sub2ind(size(RefPick),xvalBin,xvalIx));
                    xvalOffset = offsetArray(xvalIx);
                    
                    % Compute Reflected Arrival Velocity and Intercept Time
                    % IRLS Scheme
                    [xVref{rh,jj,ii}(kk,1), xToRef{rh,jj,ii}(kk,1), xDepth{rh,jj,ii}(kk,1)] ...
                        = VrmsIrls(xvalOffset,xvalRefPick);
                    % OLS Scheme
%                 [xVref{ii,jj}(kk,1), xToRef{ii,jj}(kk,1), xDepth{ii,jj}(kk,1)] ...
%                     = Vrms(xvalOffset,xvalRefPick);
                    end
                end
                if isSingleShot
                    for kk = 1:250 % 250 Random Draws
                        nCut = randsample([0,1,2],1);
                        cutChan = randsample(liveChan,nCut);
                        xvalChan = liveChan;
                        cutIx = find(ismember(liveChan,cutChan));
                        xvalChan(cutIx) = [];
                        xvalIx = find(ismember(liveChan,xvalChan));
                        
                        xvalRefPick = RefPick(xvalIx);
                        xvalOffset = offsetArray(xvalIx);                    
                    % Compute Reflected Arrival Velocity and Intercept Time
                    % OLS Scheme
                    % [xVref{ii,jj}(kk,1), xToRef{ii,jj}(kk,1), xDepth{ii,jj}(kk,1)] ...
                    %   = Vrms(xvalOffset,xvalRefPick);
                    % IRLS Scheme
                    [xVref{rh,jj,ii}(kk,1), xToRef{rh,jj,ii}(kk,1), xDepth{rh,jj,ii}(kk,1)] ...
                        = VrmsIrls(xvalOffset,xvalRefPick);    
                    end
                end

            % Create Bootstrapped Population of Surface Density
            xRhoRef{rh,jj,ii} = DryCrim([xVref{rh,jj,ii}]);
            
            % Error Analysis
            realIx = find(real([xToRef{rh,jj,ii}]));
            VoRef{rh,jj,ii} = mean([xVref{rh,jj,ii}(realIx)]);
            VoVarRef{rh,jj,ii} = var([xVref{rh,jj,ii}(realIx)]);
            toRef{rh,jj,ii} = mean([xToRef{rh,jj,ii}(realIx)]);
            toVarRef{rh,jj,ii} = var([xToRef{rh,jj,ii}(realIx)]);
            HorizonDepth{rh,jj,ii} = mean([xDepth{rh,jj,ii}(realIx)]);
            RhoRef{rh,jj,ii} = mean([xRhoRef{rh,jj,ii}(realIx)]);
            CovZP = cov([xDepth{rh,jj,ii}(realIx)],[xRhoRef{rh,jj,ii}(realIx)]);
            HorizonDepthVar{rh,jj,ii} = CovZP(1,1);
            RhoRefVar{rh,jj,ii} = CovZP(2,2);
            CovDepthRho{rh,jj,ii} = CovZP(1,2);
            
        end
        
        % Perform Dix Inversion of RMS velocities
        parfor (jj = looper, nWorkers)
%         for jj = looper
            if rh > 1
                % Real Solution Constraint
                realIxI = find(real([xToRef{rh-1,jj,ii}]));
                realIxJ = find(real([xToRef{rh,jj,ii}]));
                
                [xVint{rh,jj,ii}, xHint{rh,jj,ii}] ...
                    = DixHVA([xVref{rh-1,jj,ii}(realIxI)],[xVref{rh,jj,ii}(realIxJ)],...
                    [xToRef{rh-1,jj,ii}(realIxI)],[xToRef{rh,jj,ii}(realIxJ)]);
                xRhoInt{rh,jj,ii} = DryCrim(xVint{rh,jj,ii});
                VoInt{rh,jj,ii} = mean([xVint{rh,jj,ii}]);
                VoIntVar{rh,jj,ii} = var([xVint{rh,jj,ii}]);
                HoInt{rh,jj,ii} = mean([xHint{rh,jj,ii}]);
                HoIntVar{rh,jj,ii} = var([xHint{rh,jj,ii}]);
                RhoInt{rh,jj,ii} = mean([xRhoInt{rh,jj,ii}]);
                RhoIntVar{rh,jj,ii} = var([xRhoInt{rh,jj,ii}]);
                CovHP = cov([xHint{rh,jj,ii}],[xRhoInt{rh,jj,ii}]);
                CovThicknessRho{rh,jj,ii} = CovHP(1,2);
                
            else
                realIx = find(real([xToRef{rh,jj,ii}]));
                xVint{rh,jj,ii} = xVref{rh,jj,ii}(realIx);
                xHint{rh,jj,ii} = xDepth{rh,jj,ii}(realIx);
                xRhoInt{rh,jj,ii} = xRhoRef{rh,jj,ii}(realIx);
                HoInt{rh,jj,ii} = HorizonDepth{rh,jj,ii};
                HoIntVar{rh,jj,ii} = HorizonDepthVar{rh,jj,ii};
                VoInt{rh,jj,ii} = VoRef{rh,jj,ii};
                VoIntVar{rh,jj,ii} = VoVarRef{rh,jj,ii};
                RhoInt{rh,jj,ii} = RhoRef{rh,jj,ii};
                RhoIntVar{rh,jj,ii} = RhoRefVar{rh,jj,ii};
                CovHP = cov([xHint{rh,jj,ii} ],[xRhoInt{rh,jj,ii}]);
                CovThicknessRho{rh,jj,ii} = CovHP(1,2);
            end
            
        end

                    
        % Concatenate Reflection Travel Time
        ReflectionTo{rh,ii} = [toRef{rh,:,ii}];
        ReflectionToVar{rh,ii} = [toVarRef{rh,:,ii}];
        % Concatenate Reflection Velocity
        ReflectionVelocity{rh,ii} = [VoRef{rh,:,ii}];
        ReflectionVelocityVar{rh,ii} = [VoVarRef{rh,:,ii}];
        IntervalVelocity{rh,ii} = [VoInt{rh,:,ii}];
        IntervalVelocityVar{rh,ii} = [VoIntVar{rh,:,ii}];
        % Concatenate Reflector Depth
        ReflectionDepth{rh,ii} = [HorizonDepth{rh,:,ii}];
        ReflectionDepthVar{rh,ii} = [HorizonDepthVar{rh,:,ii}];
        %Concatenate Layer Thicnkess
        IntervalThickness{rh,ii} = [HoInt{rh,:,ii}];
        IntervalThicknessVar{rh,ii} = [HoIntVar{rh,:,ii}];
        % Concatenate Reflection Density
        ReflectionDensity{rh,ii} = [RhoRef{rh,:,ii}];
        ReflectionDensityVar{rh,ii} = [RhoRefVar{rh,:,ii}];
        % Concatenate Interval Density
        IntervalDensity{rh,ii} = [RhoInt{rh,:,ii}];
        IntervalDensityVar{rh,ii} = [RhoIntVar{rh,:,ii}];
        % Concatenate Covariance
        CovarianceDepthDensity{rh,ii} = [CovDepthRho{rh,:,ii}];
        CovarianceThicknessDensity{rh,ii} = [CovThicknessRho{rh,:,ii}];
        
        % Estimte SWE
        SWE{rh,ii} = ReflectionDepth{rh,ii}.*ReflectionDensity{rh,ii};
        SWEint{rh,ii} = IntervalThickness{rh,ii}.*IntervalDensity{rh,ii};
        % Error Propagation Equation
        SWEvar{rh,ii} = ReflectionDepthVar{rh,ii}.*ReflectionDensity{rh,ii}.^2 ...
            + ReflectionDensityVar{rh,ii}.*ReflectionDepth{rh,ii}.^2 ...
            + 2.*ReflectionDepth{rh,ii}.*ReflectionDensity{rh,ii}.*CovarianceDepthDensity{rh,ii};
        SWEintVar{rh,ii} = IntervalThicknessVar{rh,ii}.*IntervalDensity{rh,ii}.^2 ...
            + IntervalDensityVar{rh,ii}.*IntervalThickness{rh,ii}.^2 ...
            + 2.*IntervalThickness{rh,ii}.*IntervalDensity{rh,ii}.*CovarianceThicknessDensity{rh,ii};
        
        % Smooth Esitmates
        smoothR = 251;
        TWT{rh,ii} = nonParametricSmooth( 1:length(ReflectionTo{rh,ii}),...
            ReflectionTo{rh,ii},1:length(ReflectionTo{rh,ii}),smoothR);
        TWTvar{rh,ii} = nonParametricSmooth( 1:length(ReflectionToVar{rh,ii}),...
            ReflectionToVar{rh,ii},1:length(ReflectionToVar{rh,ii}),smoothR);
        SnowWaterEqv{rh,ii} = nonParametricSmooth( 1:length(SWE{rh,ii}),SWE{rh,ii},...
            1:length(SWE{rh,ii}),smoothR);
        SnowWaterEqvVar{rh,ii} = nonParametricSmooth( 1:length(SWEvar{rh,ii}),...
            SWEvar{rh,ii},1:length(SWEvar{rh,ii}),smoothR);
        Density{rh,ii} = nonParametricSmooth( 1:length(ReflectionDensity{rh,ii}),...
            ReflectionDensity{rh,ii},1:length(ReflectionDensity{rh,ii}),smoothR);
        DensityVar{rh,ii} = nonParametricSmooth( 1:length(ReflectionDensityVar{rh,ii}),...
            ReflectionDensityVar{rh,ii},1:length(ReflectionDensityVar{rh,ii}),smoothR);
        Depth{rh,ii} = nonParametricSmooth( 1:length(ReflectionDepth{rh,ii}),...
            ReflectionDepth{rh,ii},1:length(ReflectionDepth{rh,ii}),smoothR);
        DepthVar{rh,ii} = nonParametricSmooth( 1:length(ReflectionDepthVar{rh,ii}),...
            ReflectionDepthVar{rh,ii},1:length(ReflectionDepthVar{rh,ii}),smoothR);
        % Smooth Inteval Estimates
        LayerSnowWaterEqv{rh,ii} = nonParametricSmooth( 1:length(SWEint{rh,ii}),SWEint{rh,ii},...
            1:length(SWEint{rh,ii}),smoothR);
        LayerSnowWaterEqvVar{rh,ii} = nonParametricSmooth( 1:length(SWEintVar{rh,ii}),...
            SWEintVar{rh,ii},1:length(SWEintVar{rh,ii}),smoothR);
        LayerDensity{rh,ii} = nonParametricSmooth( 1:length(IntervalDensity{rh,ii}),...
            IntervalDensity{rh,ii},1:length(IntervalDensity{rh,ii}),smoothR);
        LayerDensityVar{rh,ii} = nonParametricSmooth( 1:length(IntervalDensityVar{rh,ii}),...
            IntervalDensityVar{rh,ii},1:length(IntervalDensityVar{rh,ii}),smoothR);
        LayerThickness{rh,ii} = nonParametricSmooth( 1:length(IntervalThickness{rh,ii}),...
            IntervalThickness{rh,ii},1:length(IntervalThickness{rh,ii}),smoothR);
        LayerThicknessVar{rh,ii} = nonParametricSmooth( 1:length(IntervalThicknessVar{rh,ii}),...
            IntervalThicknessVar{rh,ii},1:length(IntervalThicknessVar{rh,ii}),smoothR);
        end
    end
    
    fprintf('Horizon Velocity Analysis Done \n')
    toc
    display(' ')

else    % If No Time Correction is Used ~ Not Recommended
    for ii = 1:nFiles
        % Record Time
        time{ii} = [0:dt:(size(Radar{ii,1},1)-(1+padding(ii)))*dt];
    end
end

%% Save HVA Output
isSaveHVA = 0;
if isSaveHVA
    % Save Rough Results
    HVAfilename = '6-2-16-Core7-Spur-W-HVA.mat';
    save(HVAfilename,'DirectTo','DirectToVar','deltaT',...
    'DirectVelocity','DirectVelocityVar','DirectDepth','DirectDepthVar',...
    'DirectDensity','DirectDensityVar','CovarianceDepthDensityDirect',...
    'DirectSWE','DirectSWEvar','ReflectionTo','ReflectionToVar',...
    'ReflectionVelocity','ReflectionVelocityVar','IntervalVelocity',...
    'IntervalVelocityVar','ReflectionDepth','ReflectionDepthVar',...
    'IntervalThickness','IntervalThicknessVar','ReflectionDensity',...
    'ReflectionDensityVar','IntervalDensity','IntervalDensityVar',...
    'CovarianceDepthDensity','CovarianceThicknessDensity','SWE','SWEint',...
    'SWEvar','SWEintVar','-v7.3');

    % Save Smooth Results
    HVAsmoothFilename = '6-2-16-Core7-Spur-W-HVA.mat';
    save(HVAsmoothFilename,'TWT','TWTvar','SnowWaterEqv',...
        'SnowWaterEqvVar','Density','DensityVar','Depth','DepthVar','LayerSnowWaterEqv',...
        'LayerSnowWaterEqvVar','LayerDensity','LayerDensityVar','LayerThickness',...
        'LayerThicknessVar','-v7.3');
end
%% Create Figures for Snow Surface Data
if isSWEDISH
    % plot Direct Wave Data
    for ii = 1:nFiles
        figure();
        subplot(3,1,1)
        for dh = 2:nDirectHorizon
            shadedErrorBarT8([],dhSnowWaterEqv{dh,ii},...
                sqrt(dhSnowWaterEqvVar{dh,ii}),1,{'Color',[0.5,0,0],'linewidth',1.5});
            hold on;
        end
        freezeColors
        axis tight
        grid on
        title('Accumulation [m w.e.]')
        subplot(3,1,2)
        for dh = 2:nDirectHorizon
            shadedErrorBarT8([],dhDensity{dh,ii}.*1000,...
                sqrt(dhDensityVar{dh,ii}).*1000,1,{'Color',[0,0,.5],'linewidth',1.5});
            hold on;
        end
        freezeColors
        axis tight
        grid on
        title('Density [kgm^{-3}]')
        subplot(3,1,3)
        for dh = 2:nDirectHorizon
            shadedErrorBarT8([],dhDepth{dh,ii},...
                sqrt(dhDepthVar{dh,ii}),1,{'Color',[0,0,0],'linewidth',1.5});
            hold on;
        end
        freezeColors
        axis tight
        grid on
        title('Depth [m]')
        set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,...
            'FontWeight','Bold', 'LineWidth', 1);
    end
    % Plot Reflection Data
    for ii = 1:nFiles
        figure();
        subplot(3,1,1)
        for rh = 1:nReflectionHorizon
            shadedErrorBarT8([],SnowWaterEqv{rh,ii},...
                sqrt(SnowWaterEqvVar{rh,ii}),1,{'Color',[0.5,0,0],'linewidth',1.5});
            hold on;
        end
        freezeColors
        axis tight
        grid on
        title('Accumulation [m w.e.]')
        subplot(3,1,2)
        for rh = 1:nReflectionHorizon
            shadedErrorBarT8([],Density{rh,ii}.*1000,...
                sqrt(DensityVar{rh,ii}).*1000,1,{'Color',[0,0,.5],'linewidth',1.5});
            hold on;
        end
        freezeColors
        axis tight
        grid on
        title('Density [kgm^{-3}]')
        subplot(3,1,3)
        for rh = 1:nReflectionHorizon
            shadedErrorBarT8([],Depth{rh,ii},...
                sqrt(DepthVar{rh,ii}),1,{'Color',[0,0,0],'linewidth',1.5});
            hold on;
        end
        freezeColors
        axis tight
        grid on
        title('Depth [m]')
        set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,...
            'FontWeight','Bold', 'LineWidth', 1);
    end
end