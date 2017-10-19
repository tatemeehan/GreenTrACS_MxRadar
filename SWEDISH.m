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

% Process Data
isTrimTWT = 0;          % Truncate Recorded Data for Near-Surface Analysis
isKill = 0;             % Kill Unanted Channels
isFXdecon = 1;          % Fx-Predicitive Deconvolution
isSWEDISH = 1;          % Perform Surface Velocity Analysis

% Export Data
isWrite = 0;

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
            trhd{ii}(1,jj:end) = trhd{ii}(1,jj:end) - 1; % Configure Trace Indicies
            trhd{ii}(2,jj:end) = trhd{ii}(2,jj:end) - dx;% Configure Distance
        end
        
        trhd{ii}(:,dupIx) = []; % Remove Skipped Traces from Trace Header
        Rad{ii}(:,dupIx) = []; % Remove Skipped Traces from Multiplexed Data
        xArray = trhd{ii}(2,:); % Define Configured Distance as xArray
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
    xcorrWindow = [175, 250, 85, 140, 230, 75, 100, 180, 25;...
                    200, 300, 150, 160, 260, 115, 150, 220, 75] + padding;
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
%% FX Predicive Deconvolution
if isFXdecon
    fprintf('Begin FX-Predictive Deconvolution \n')
    tic
    saveRadar = Radar;
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
    %% Semi-Automatic Radar Wave Picking
    fprintf('Begin PolarPicker \n')
    display(' ')
    tic
    
    pickRadar = Radar;

    % AGC Gain for PolarPicker 
    for ii = 1:nFiles
        parfor (jj = chan, nWorkers)
            pickRadar{jj,ii} = AGCgain(pickRadar{jj,ii},50,2);
        end
    end
    
    % Pick Direct Wave Arrival
    display('Pick Direct Wave')
    display(' ')
%     display('If Multiple Horizons are Chosen, Initialize with Air Wave')
%     display(' ')
    [~, DirectFBpick] = polarPicker(pickRadar);
    
    pickRadar = Radar;
    
    % Median Subtraction Filter For Reflection Analysis
    for ii = 1:nFiles
        parfor (jj = chan, nWorkers)
            pickRadar{jj,ii} = pickRadar{jj,ii}-median(pickRadar{jj,ii},2)...
                *ones(1,size(pickRadar{jj,ii},2));
        end
    end

    % AGC Gain for PolarPicker 
    for ii = 1:nFiles
        parfor (jj = chan, nWorkers)
            pickRadar{jj,ii} = AGCgain(pickRadar{jj,ii},350,2);
        end
    end
    
    % Pick Primary Reflection Arrival
    display('Pick Primary Reflection')
    display(' ')
    [~, ReflectionFBpick] = polarPicker(pickRadar);
    
    clear pickRadar
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

    fprintf('PolarPicker Done \n')
    toc
    display(' ')
        
%% Perform Horizon Velocity Analysis for Estimates of Density, Depth, & SWE
    fprintf('Begin Horizon Velocity Analysis \n')
    display(' ')
    tic
  
    % Determine Maximum File Size
    nTrace = zeros(1,nFiles);
    for ii = 1:nFiles
        nTrace(ii) = size(Radar{1,ii},2);
    end
    % Maxiumum Traces to Loop Over
    nTrace = max(nTrace);
        
    % Allocate Memory
    GatherReflectionPicks = cell(nReflectionHorizon,nFiles);
    GatherDirectPicks = cell(nDirectHorizon,nFiles);
    
    % Allocate Direct Wave LMO Velocity Boothstrapping    
    xVdir = cell(nDirectHorizon,nTrace,nFiles);xToDir = cell(nDirectHorizon,nTrace,nFiles);
    VoDir = cell(nDirectHorizon,nTrace,nFiles);VoVarDir = cell(nDirectHorizon,nTrace,nFiles);
    toDir = cell(nDirectHorizon,nTrace,nFiles);toVarDir = cell(nDirectHorizon,nTrace,nFiles);
    DirVelocity = cell(nDirectHorizon,nFiles);RhoDir = cell(nDirectHorizon,nTrace,nFiles);
    xRhoDir = cell(nDirectHorizon,nTrace,nFiles);DirDensity = cell(nDirectHorizon,nFiles);
    xDepthDir = cell(nDirectHorizon,nTrace,nFiles);DepthDir  = cell(nDirectHorizon,nTrace,nFiles);
    DepthVarDir = cell(nDirectHorizon,nTrace,nFiles);RhoDirVar  = cell(nDirectHorizon,nTrace,nFiles);
    CovDepthRhoDir = cell(nDirectHorizon,nTrace,nFiles);
    deltaT = cell(nFiles,nTrace);ResGatherDirPicks = cell(nDirectHorizon,nTrace,nFiles);
    AirTo = cell(nFiles,nTrace); xAirTo = cell(nFiles,nTrace);
    
    % Allocate Direct Wave Snow Analysis
    DirectTo = cell(nDirectHorizon,nFiles);DirectToVar = cell(nDirectHorizon,nFiles);
    DirectVelocity = cell(nDirectHorizon,nFiles); DirectVelocityVar = cell(nDirectHorizon,nFiles);
    DirectDepth = cell(nDirectHorizon,nFiles); DirectDepthVar = cell(nDirectHorizon,nFiles);
    DirectDensity = cell(nDirectHorizon,nFiles); DirectDensityVar = cell(nDirectHorizon,nFiles);
    CovarianceDepthDensityDirect = cell(nDirectHorizon,nFiles); 
    DirectSWE = cell(nDirectHorizon,nFiles); DirectSWEvar = cell(nDirectHorizon,nFiles);
    
    % Allocate for Concatenation and Plotting
    dhTWT = cell(nDirectHorizon,nFiles); dhTWTvar = cell(nDirectHorizon,nFiles);
    dhSnowWaterEqv = cell(nDirectHorizon,nFiles);dhSnowWaterEqvVar = cell(nDirectHorizon,nFiles);
    dhDensity = cell(nDirectHorizon,nFiles); dhDensityVar = cell(nDirectHorizon,nFiles); 
    dhDepth = cell(nDirectHorizon,nFiles); dhDepthVar = cell(nDirectHorizon,nFiles);
   
    % Allocate Bootstrapping of RMS Velocity
    xVref = cell(nReflectionHorizon,nTrace,nFiles);xToRef = cell(nReflectionHorizon,nTrace,nFiles);
    xDepth = cell(nReflectionHorizon,nTrace,nFiles);xRhoRef = cell(nReflectionHorizon,nTrace,nFiles);
    VoRef = cell(nReflectionHorizon,nTrace,nFiles);VoVarRef = cell(nReflectionHorizon,nTrace,nFiles);
    toRef = cell(nReflectionHorizon,nTrace,nFiles);toVarRef = cell(nReflectionHorizon,nTrace,nFiles);
    HorizonDepth = cell(nReflectionHorizon,nTrace,nFiles);RhoRef = cell(nReflectionHorizon,nTrace,nFiles);
    HorizonDepthVar = cell(nReflectionHorizon,nTrace,nFiles);RhoRefVar = cell(nReflectionHorizon,nTrace,nFiles);
    CovDepthRho = cell(nReflectionHorizon,nTrace,nFiles);
    
    % Allocation for Interval Velocities
    xVint = cell(nReflectionHorizon,nTrace,nFiles);xHint = cell(nReflectionHorizon,nTrace,nFiles);
    VoInt = cell(nReflectionHorizon,nTrace,nFiles);VoIntVar = cell(nReflectionHorizon,nTrace,nFiles);
    HoInt = cell(nReflectionHorizon,nTrace,nFiles);HoIntVar = cell(nReflectionHorizon,nTrace,nFiles);
    xRhoInt = cell(nReflectionHorizon,nTrace,nFiles); RhoInt = cell(nReflectionHorizon,nTrace,nFiles);
    RhoIntVar = cell(nReflectionHorizon,nTrace,nFiles);CovThicknessRho = cell(nReflectionHorizon,nTrace,nFiles);      
    
    % Allocate for Concatenation and Plotting    
    ReflectionVelocity = cell(nReflectionHorizon,nFiles);
    ReflectionVelocityVar = cell(nReflectionHorizon,nFiles);ReflectionDepth = cell(nReflectionHorizon,nFiles);
    ReflectionDepthVar = cell(nReflectionHorizon,nFiles);ReflectionDensity = cell(nReflectionHorizon,nFiles);
    ReflectionDensityVar = cell(nReflectionHorizon,nFiles);CovarianceDepthDensity = cell(nReflectionHorizon,nFiles);
    ReflectionTo = cell(nReflectionHorizon,nFiles); ReflectionToVar = cell(nReflectionHorizon,nFiles);
    IntervalVelocity = cell(nReflectionHorizon,nFiles);IntervalVelocityVar = cell(nReflectionHorizon,nFiles);
    IntervalThickness = cell(nReflectionHorizon,nFiles);IntervalThicknessVar = cell(nReflectionHorizon,nFiles);
    IntervalDensity = cell(nReflectionHorizon,nFiles);IntervalDensityVar = cell(nReflectionHorizon,nFiles);
    CovarianceThicknessDensity = cell(nReflectionHorizon,nFiles);
    
    % Allocation for Results
    SWE = cell(nReflectionHorizon,nFiles);SWEvar = cell(nReflectionHorizon,nFiles);
    SWEint = cell(nReflectionHorizon,nFiles);SWEintVar = cell(nReflectionHorizon,nFiles);
    TWT = cell(nReflectionHorizon,nFiles);TWTvar = cell(nReflectionHorizon,nFiles);
    SnowWaterEqv = cell(nReflectionHorizon,nFiles);SnowWaterEqvVar = cell(nReflectionHorizon,nFiles);
    Density = cell(nReflectionHorizon,nFiles); DensityVar = cell(nReflectionHorizon,nFiles);
    Depth = cell(nReflectionHorizon,nFiles); DepthVar = cell(nReflectionHorizon,nFiles);
    LayerSnowWaterEqv = cell(nReflectionHorizon,nFiles);LayerSnowWaterEqvVar = cell(nReflectionHorizon,nFiles);
    LayerDensity = cell(nReflectionHorizon,nFiles);LayerDensityVar = cell(nReflectionHorizon,nFiles);
    LayerThickness = cell(nReflectionHorizon,nFiles);LayerThicknessVar = cell(nReflectionHorizon,nFiles);    
    
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
                    xDepthDir{dh,jj,ii} = xVdir{dh,jj,ii}.*xToDir{dh,jj,ii};
                    % Estimate Direct Wave Density
                    xRhoDir{dh,jj,ii} = DryCrim([xVdir{dh,jj,ii}]);
                    % Error Analysis
                    VoDir{dh,jj,ii} = mean([xVdir{dh,jj,ii}]);
                    VoVarDir{dh,jj,ii} = var([xVdir{dh,jj,ii}]);
                    toDir{dh,jj,ii} = mean([xToDir{dh,jj,ii}]);
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
%% Old Version 9 - 5 - 17    
%     %% Semi-Automatic Radar Wave Picking
%     fprintf('Begin PolarPicker \n')
%     display(' ')
%     tic
%     
%     pickRadar = Radar;
% 
%     % AGC Gain for PolarPicker 
%     for ii = 1:nFiles
%         parfor (jj = chan, nWorkers)
%             pickRadar{jj,ii} = AGCgain(pickRadar{jj,ii},50,2);
%         end
%     end
%     
%     % Pick AirWave Arrival
%     display('Pick Air Wave')
%     display(' ')
%     [~, AirFBpick] = polarPicker(pickRadar);
%     
%     % Pick Primary Reflection Arrival
%     display('Pick Primary Reflection')
%     display(' ')
%     [~, ReflectionFBpick] = polarPicker(pickRadar);
%     
%     clear pickRadar
%     
%     % Determine Number of Horizons
%     nHorizon = size(ReflectionFBpick,2);
%   
% %     % Apply Systematic Correction to First Break Picks
% %         time = cell(1,nFiles);
% %         for ii = 1:nFiles
% %             for jj = 1:nChan
% %                 ReflectionFBpick{jj,ii} = (ReflectionFBpick{jj,ii} - chanShift(jj));
% %                 AirFBpick{jj,ii} = (AirFBpick{jj,ii} - chanShift(jj));
% %             end
% %         end
%         
% % Convert Samples to ns
% for ii = 1:nFiles
%     time{ii} = [0:dt:(size(Radar{ii,1},1)-(1+padding(ii)))*dt];
%     for jj = chan
%         AirFBpick{jj,ii} = AirFBpick{jj,ii}.*dt;
%         for hh = 1:nHorizon
%             ReflectionFBpick{jj,hh,ii} = ReflectionFBpick{jj,hh,ii}.*dt;
%         end
%     end
% end
% 
%     fprintf('PolarPicker Done \n')
%     toc
%     display(' ')
%         
% %% Perform Horizon Velocity Analysis for Estimates of Density, Depth, & SWE
% % First Horizon has Bad to and depth Complex. Perhaps Bad Picks
%     fprintf('Begin Surface Velocity Analysis \n')
%     display(' ')
%     tic
%     
%     % Determin Maximum File Size
%     for ii = 1:nFiles
%         nTrace(ii) = length(Radar{1,ii});
%     end
%     nTrace = max(nTrace);
%         
%     % Allocate Memory
%     GatherReflectionPicks = cell(nHorizon,nFiles);
%     GatherAirPicks = cell(1,nFiles);
%     
%     xVair = cell(nFiles,nTrace);xToAir = cell(nFiles,nTrace);
%     VoAir = cell(nFiles,nTrace);VoVarAir = cell(nFiles,nTrace);
%     toAir = cell(nFiles,nTrace);toVarAir = cell(nFiles,nTrace);
%     AirVelocity = cell(nFiles,1);AirDensity = cell(nFiles,1); 
%     deltaTair = cell(nFiles,nTrace);ResGatherAirPicks = cell(nFiles,nTrace);
%    
%     xVref = cell(nHorizon,nTrace,nFiles);xToRef = cell(nHorizon,nTrace,nFiles);
%     xDepth = cell(nHorizon,nTrace,nFiles);xRhoRef = cell(nHorizon,nTrace,nFiles);
%     VoRef = cell(nHorizon,nTrace,nFiles);VoVarRef = cell(nHorizon,nTrace,nFiles);
%     toRef = cell(nHorizon,nTrace,nFiles);toVarRef = cell(nHorizon,nTrace,nFiles);
%     HorizonDepth = cell(nHorizon,nTrace,nFiles);RhoRef = cell(nHorizon,nTrace,nFiles);
%     HorizonDepthVar = cell(nHorizon,nTrace,nFiles);RhoRefVar = cell(nHorizon,nTrace,nFiles);
%     CovDepthRho = cell(nHorizon,nTrace,nFiles);ReflectionVelocity = cell(nHorizon,nFiles);
%     ReflectionVelocityVar = cell(nHorizon,nFiles);ReflectionDepth = cell(nHorizon,nFiles);
%     ReflectionDepthVar = cell(nHorizon,nFiles);ReflectionDensity = cell(nHorizon,nFiles);
%     ReflectionDensityVar = cell(nHorizon,nFiles);CovarianceDepthDensity = cell(nHorizon,nFiles);
%     ReflectionTo = cell(nHorizon,nFiles); ReflectionToVar = cell(nHorizon,nFiles);
%     SWE = cell(nHorizon,nFiles);SWEvar = cell(nHorizon,nFiles);
%     TWT = cell(nHorizon,nFiles);TWTvar = cell(nHorizon,nFiles);
%     SnowWaterEqv = cell(nHorizon,nFiles);SnowWaterEqvVar = cell(nHorizon,nFiles);
%     Density = cell(nHorizon,nFiles); DensityVar = cell(nHorizon,nFiles);
%     Depth = cell(nHorizon,nFiles); DepthVar = cell(nHorizon,nFiles);
%     
% % Extract AirWave Picks and Perform Residual Subtraction Velocity Analysis
% %     parfor (ii = 1:nFiles, nWorkers - (nWorkers-nFiles)) 
%     for ii = 1:nFiles
%         looper = 1:length(Radar{ii});
%         % Concatenate Air Wave Picks
%         GatherAirPicks{ii} = cat(2,AirFBpick{:,ii});
%         Air = GatherAirPicks{ii};
%         % Air Wave Velocity Analysis      
%         parfor (jj = looper, nWorkers)
% %         for jj = 1:length(Radar{ii})
%             % Extract Picks
%             AirPick = Air(jj,:); 
%             
%             for kk = 1:100 % 100 Random Draws
%                 % Cross-Validation for Surface Velocity Estimation 7-28-17
%                 nCut = randsample([0,1,2],1);
%                 cutChan = randsample(liveChan,nCut);
%                 xvalChan = liveChan;
%                 cutIx = find(ismember(liveChan,cutChan));
%                 xvalChan(cutIx) = [];
%                 xvalIx = find(ismember(xvalChan,liveChan));
%                 
%                 xvalAirPick = AirPick(xvalIx);
%                 xvalOffset = offsetArray(xvalIx);
%                 
%                 % Compute Air Wave Arrival Velocity and Intercept Time
%                 % OLS Scheme
% %                 [xVair{ii,jj}(kk,1), ~] = DirectWave(xvalOffset,xvalAirPick);
%                 % IRLS Scheme
%                 [xVair{ii,jj}(kk,1), xToAir{ii,jj}(kk,1)] = DirectWaveIrls(xvalOffset,xvalAirPick);                
%             end
% 
%             % Residual Error Analysis 
%             VoAir{ii,jj} = mean([xVair{ii,jj}]);
% 
%             % Find Velocity Residual for Additional Trace Shifts & Calulate Residual Time
%             deltaTair{ii,jj} = (offsetArray./VoAir{ii,jj}) - (offsetArray./0.299);
% 
%             % Apply Residual Time Shifts
%             ResGatherAirPicks{ii,jj} = Air(jj,:) - deltaTair{ii,jj};
% 
%             % Re-Extract Airwave Picks after Residual Shifts
%             AirPick = ResGatherAirPicks{ii,jj}; 
%             
%             for kk = 1:100 % 100 Random Draws
%                 % Cross-Validation for Surface Velocity Estimation 7-28-17
%                 nCut = randsample([0,1,2],1);
%                 cutChan = randsample(liveChan,nCut);
%                 xvalChan = liveChan;
%                 cutIx = find(ismember(liveChan,cutChan));
%                 xvalChan(cutIx) = [];                
%                 xvalIx = find(ismember(xvalChan,liveChan));
%                 
%                 xvalAirPick = AirPick(xvalIx);
%                 xvalOffset = offsetArray(xvalIx);
%                 
%                 % Compute Air Wave Arrival Velocity and Intercept Time
%                 % OLS Scheme
% %                 [xVair{ii,jj}(kk,1), xToAir{ii,jj}(kk,1)] = DirectWave(xvalOffset,xvalAirPick);
%                 % IRLS Scheme
%                 [xVair{ii,jj}(kk,1), xToAir{ii,jj}(kk,1)] = DirectWaveIrls(xvalOffset,xvalAirPick);                
%             end
%             % Error Analysis 
%             VoAir{ii,jj} = mean([xVair{ii,jj}]);
%             VoVarAir{ii,jj} = var([xVair{ii,jj}]);
%             toAir{ii,jj} = mean([xToAir{ii,jj}]);
%             toVarAir{ii,jj} = var([xToAir{ii,jj}]);
%         end
%         % Concatenate Air Velocity
%         AirVelocity{ii} = [VoAir{ii,:}];
%         % Estimate Air Density
%         AirDensity{ii} = DryCrim(AirVelocity{ii});
%     end
%             
% 
% % Primary Reflection Velocity Analysis               
%     for ii = 1:nFiles
%         looper = 1:length(Radar{ii});
%         for hh = 1:nHorizon
%         % Concatenate Primary Reflection Picks for Horizon hh
%         GatherReflectionPicks{hh,ii} = cat(2,ReflectionFBpick{:,hh,ii});
%         Reflection = GatherReflectionPicks{hh,ii};
% %         for jj = 1:length(Radar{ii})
%         parfor (jj = looper, nWorkers)
%             % Cross-Validation for Reflection Velocity Estimation
%             % Multiple Shot Gathers in Population
%             isManyShots = 0;
%             if isManyShots
%                 range = sqrt(([1:length(Radar{ii})]-jj).^2); % Compute Distance
%                 getIx = find(range<=5); % Find Nearby Picks
%                 RefPick = GatherReflectionPicks{ii}(getIx,:) - vertcat(deltaTair{ii,getIx}); % Residual
%                 RefPick = RefPick - vertcat(toAir{ii,getIx})*ones(1,nChan);
%                 stdRefPick = std(RefPick,0,1);
%             end
%             % Single Shot Gather in Population
%             isSingleShot = 1;
%             if isSingleShot
%                 RefPick = Reflection(jj,:) - deltaTair{ii,jj}; % Residual
%                 RefPick = RefPick - toAir{ii,jj};
%             end
% 
%             for kk = 1:100 % 100 Random Draws
%                 % Cross-Validation for Surface Velocity Estimation 7-28-17
%                 nCut = randsample([0,1,2],1);
%                 cutChan = randsample(liveChan,nCut);
%                 xvalChan = liveChan;
%                 cutIx = find(ismember(liveChan,cutChan));
%                 xvalChan(cutIx) = [];                
%                 xvalIx = find(ismember(xvalChan,liveChan));
%                 
%                 xvalRefPick = RefPick(xvalIx);
%                 xvalOffset = offsetArray(xvalIx);
%                 
%                 % Compute Reflected Arrival Velocity and Intercept Time
%                 % OLS Scheme
% %                 [xVref{ii,jj}(kk,1), xToRef{ii,jj}(kk,1), xDepth{ii,jj}(kk,1)] ...
% %                     = Vrms(xvalOffset,xvalRefPick);
%                 % IRLS Scheme
% %                 [xVref{hh,jj,ii}(kk,1), xToRef{hh,jj,ii}(kk,1), xDepth{hh,jj,ii}(kk,1)] ...
% %                     = VrmsIrls(xvalOffset,xvalRefPick);
%                 % Linear Arrival Approach
%                 [xVref{hh,jj,ii}(kk,1), xToRef{hh,jj,ii}(kk,1)] ...
%                     = DirectWaveIrls(xvalOffset,xvalRefPick);
% %                 xToRef{hh,jj,ii} = abs(xToRef{hh,jj,ii}(kk,1));
%                 xDepth{hh,jj,ii}(kk,1) = xVref{hh,jj,ii}(kk,1).*xToRef{hh,jj,ii}(kk,1);
%             end
% %             xToRef{hh,jj,ii} = abs(xToRef{hh,jj,ii});
% %             xDepth{hh,jj,ii} = abs(xDepth{hh,jj,ii});
%             xDepth{hh,jj,ii} = xVref{hh,jj,ii}.*xToRef{hh,jj,ii};
%             % Create Random Sample Population of Surface Density
%             xRhoRef{hh,jj,ii} = DryCrim([xVref{hh,jj,ii}]);
%             
%             % Error Analysis
%             realIx = find(real([xToRef{hh,jj,ii}]));
%             VoRef{hh,jj,ii} = mean([xVref{hh,jj,ii}(realIx)]);
%             VoVarRef{hh,jj,ii} = var([xVref{hh,jj,ii}(realIx)]);
%             toRef{hh,jj,ii} = mean([xToRef{hh,jj,ii}(realIx)]);
%             toVarRef{hh,jj,ii} = var([xToRef{hh,jj,ii}(realIx)]);
%             HorizonDepth{hh,jj,ii} = mean([xDepth{hh,jj,ii}(realIx)]);
%             RhoRef{hh,jj,ii} = mean([xRhoRef{hh,jj,ii}(realIx)]);
%             CovZP = cov([xDepth{hh,jj,ii}(realIx)],[xRhoRef{hh,jj,ii}(realIx)]);
%             HorizonDepthVar{hh,jj,ii} = CovZP(1,1);
%             RhoRefVar{hh,jj,ii} = CovZP(2,2);
%             CovDepthRho{hh,jj,ii} = CovZP(1,2);
%             
%         end
%         % Concatenate Reflection Travel Time
%         ReflectionTo{hh,ii} = [toRef{hh,:,ii}];
%         ReflectionToVar{hh,ii} = [toVarRef{hh,:,ii}];
%         % Concatenate Reflection Velocity
%         ReflectionVelocity{hh,ii} = [VoRef{hh,:,ii}];
%         ReflectionVelocityVar{hh,ii} = [VoVarRef{hh,:,ii}];
%         % Concatenate Reflector Depth
%         ReflectionDepth{hh,ii} = [HorizonDepth{hh,:,ii}];
%         ReflectionDepthVar{hh,ii} = [HorizonDepthVar{hh,:,ii}];
%         % Concatenate Reflection Density & Convert Units
%         ReflectionDensity{hh,ii} = [RhoRef{hh,:,ii}];
%         ReflectionDensityVar{hh,ii} = [RhoRefVar{hh,:,ii}];
%         % Concatenate Covariance
%         CovarianceDepthDensity{hh,ii} = [CovDepthRho{hh,:,ii}];
%         
%         % Estimte SWE
%         SWE{hh,ii} = ReflectionDepth{hh,ii}.*ReflectionDensity{hh,ii};
%         % Error Propagation Equation
%         SWEvar{hh,ii} = ReflectionDepthVar{hh,ii}.*ReflectionDensity{hh,ii}.^2 ...
%             + ReflectionDensityVar{hh,ii}.*ReflectionDepth{hh,ii}.^2 ...
%             + 2.*ReflectionDepth{hh,ii}.*ReflectionDensity{hh,ii}.*CovarianceDepthDensity{hh,ii};
%         
%         % Smooth Esitmates
%         smoothR = 251;
%         TWT{hh,ii} = nonParametricSmooth( 1:length(ReflectionTo{hh,ii}),...
%             ReflectionTo{hh,ii},1:length(ReflectionTo{hh,ii}),smoothR);
%         TWTvar{hh,ii} = nonParametricSmooth( 1:length(ReflectionToVar{hh,ii}),...
%             ReflectionToVar{hh,ii},1:length(ReflectionToVar{hh,ii}),smoothR);
%         SnowWaterEqv{hh,ii} = nonParametricSmooth( 1:length(SWE{hh,ii}),SWE{hh,ii},...
%             1:length(SWE{hh,ii}),smoothR);
%         SnowWaterEqvVar{hh,ii} = nonParametricSmooth( 1:length(SWEvar{hh,ii}),...
%             SWEvar{hh,ii},1:length(SWEvar{hh,ii}),smoothR);
%         Density{hh,ii} = nonParametricSmooth( 1:length(ReflectionDensity{hh,ii}),...
%             ReflectionDensity{hh,ii},1:length(ReflectionDensity{hh,ii}),smoothR);
%         DensityVar{hh,ii} = nonParametricSmooth( 1:length(ReflectionDensityVar{hh,ii}),...
%             ReflectionDensityVar{hh,ii},1:length(ReflectionDensityVar{hh,ii}),smoothR);
%         Depth{hh,ii} = nonParametricSmooth( 1:length(ReflectionDepth{hh,ii}),...
%             ReflectionDepth{hh,ii},1:length(ReflectionDepth{hh,ii}),smoothR);
%         DepthVar{hh,ii} = nonParametricSmooth( 1:length(ReflectionDepthVar{hh,ii}),...
%             ReflectionDepthVar{hh,ii},1:length(ReflectionDepthVar{hh,ii}),smoothR);
%         end
%     end
%     
%     fprintf('Surface Velocity Analysis Done \n')
%     toc
%     display(' ')
% 
% else    % If No Time Correction is Used ~ Not Recommended
%     for ii = 1:nFiles
%         % Record Time
%         time{ii} = [0:dt:(size(Radar{ii,1},1)-(1+padding(ii)))*dt];
%     end
% end

%% Old Verion
% %% Perform Surface Velocity Analysis for Estimates of Density, Depth, & SWE
% 
%     fprintf('Begin Surface Velocity Analysis \n')
%     display(' ')
%     tic
%     
%     % Determin Maximum File Size
%     for ii = 1:nFiles
%         nTrace(ii) = length(Radar{1,ii});
%     end
%     nTrace = max(nTrace);
%     
%     % Allocate Memory
%     GatherReflectionPicks = cell(1,nFiles);
%     GatherAirPicks = cell(1,nFiles);
%     
%     xVair = cell(nFiles,nTrace);xToAir = cell(nFiles,nTrace);
%     VoAir = cell(nFiles,nTrace);VoVarAir = cell(nFiles,nTrace);
%     toAir = cell(nFiles,nTrace);toVarAir = cell(nFiles,nTrace);
%    
%     xVref = cell(nFiles,nTrace);xToRef = cell(nFiles,nTrace);xDepth = cell(nFiles,nTrace);
%     xRhoRef = cell(nFiles,nTrace);VoRef = cell(nFiles,nTrace);VoVarRef = cell(nFiles,nTrace);
%     toRef = cell(nFiles,nTrace);toVarRef = cell(nFiles,nTrace);Depth = cell(nFiles,nTrace);
%     RhoRef = cell(nFiles,nTrace);DepthVar = cell(nFiles,nTrace);RhoRefVar = cell(nFiles,nTrace);
%     CovDepthRho = cell(nFiles,nTrace);ReflectionVelocity = cell(1,nFiles);
%     ReflectionVelocityVar = cell(1,nFiles);ReflectionDepth = cell(1,nFiles);
%     ReflectionDepthVar = cell(1,nFiles);ReflectionDensity = cell(1,nFiles);
%     ReflectionDensityVar = cell(1,nFiles);CovarianceDepthDensity = cell(1,nFiles);
%     ReflectionTo = cell(1,nFiles); ReflectionToVar = cell(1,nFiles);
%     SWE = cell(1,nFiles);SWEvar = cell(1,nFiles);AirVelocity = cell(nFiles,nTrace);
%     AirDensity = cell(nFiles,nTrace); deltaTair = cell(nFiles,nTrace);
%     ResGatherAirPicks = cell(nFiles,nTrace);
% 
% % Extract AirWave Picks and Perform Residual Subtraction Velocity Analysis
% %     parfor (ii = 1:nFiles, nWorkers - (nWorkers-nFiles)) 
%     for ii = 1:nFiles
%         looper = 1:length(Radar{ii});
%         % Concatenate Air Wave Picks
%         GatherAirPicks{ii} = cat(2,AirFBpick{:,ii});
%         Air = GatherAirPicks{ii};
%         % Air Wave Velocity Analysis      
%         parfor (jj = looper, nWorkers)
% %         for jj = 1:length(Radar{ii})
%             % Extract Picks
%             AirPick = Air(jj,:); 
%             
%             for kk = 1:100 % 100 Random Draws
%                 % Cross-Validation for Surface Velocity Estimation 7-28-17
%                 nCut = randsample([0,1,2],1);
%                 cutChan = randsample(liveChan,nCut);
%                 xvalChan = liveChan;
%                 cutIx = find(ismember(liveChan,cutChan));
%                 xvalChan(cutIx) = [];
%                 xvalIx = find(ismember(xvalChan,liveChan));
%                 
%                 xvalAirPick = AirPick(xvalIx);
%                 xvalOffset = offsetArray(xvalIx);
%                 
%                 % Compute Air Wave Arrival Velocity and Intercept Time
%                 % OLS Scheme
% %                 [xVair{ii,jj}(kk,1), ~] = DirectWave(xvalOffset,xvalAirPick);
%                 % IRLS Scheme
%                 [xVair{ii,jj}(kk,1), xToAir{ii,jj}(kk,1)] = DirectWaveIrls(xvalOffset,xvalAirPick);                
%             end
% 
%             % Residual Error Analysis 
%             VoAir{ii,jj} = mean([xVair{ii,jj}]);
% 
%             % Find Velocity Residual for Additional Trace Shifts & Calulate Residual Time
%             deltaTair{ii,jj} = (offsetArray./VoAir{ii,jj}) - (offsetArray./0.299);
% 
%             % Apply Residual Time Shifts
%             ResGatherAirPicks{ii,jj} = Air(jj,:) - deltaTair{ii,jj};
% 
%             % Re-Extract Airwave Picks after Residual Shifts
%             AirPick = ResGatherAirPicks{ii,jj}; 
%             
%             for kk = 1:100 % 100 Random Draws
%                 % Cross-Validation for Surface Velocity Estimation 7-28-17
%                 nCut = randsample([0,1,2],1);
%                 cutChan = randsample(liveChan,nCut);
%                 xvalChan = liveChan;
%                 cutIx = find(ismember(liveChan,cutChan));
%                 xvalChan(cutIx) = [];                
%                 xvalIx = find(ismember(xvalChan,liveChan));
%                 
%                 xvalAirPick = AirPick(xvalIx);
%                 xvalOffset = offsetArray(xvalIx);
%                 
%                 % Compute Air Wave Arrival Velocity and Intercept Time
%                 % OLS Scheme
% %                 [xVair{ii,jj}(kk,1), xToAir{ii,jj}(kk,1)] = DirectWave(xvalOffset,xvalAirPick);
%                 % IRLS Scheme
%                 [xVair{ii,jj}(kk,1), xToAir{ii,jj}(kk,1)] = DirectWaveIrls(xvalOffset,xvalAirPick);                
%             end
%             % Error Analysis 
%             VoAir{ii,jj} = mean([xVair{ii,jj}]);
%             VoVarAir{ii,jj} = var([xVair{ii,jj}]);
%             toAir{ii,jj} = mean([xToAir{ii,jj}]);
%             toVarAir{ii,jj} = var([xToAir{ii,jj}]);
%         end
%         % Concatenate Air Velocity
%         AirVelocity{ii} = [VoAir{ii,:}];
%         % Estimate Air Density
%         AirDensity{ii} = DryCrim(AirVelocity{ii});
%     end
%             
% 
% % Primary Reflection Velocity Analysis               
%     for ii = 1:nFiles
%         looper = 1:length(Radar{ii});
%         % Concatenate Primary Reflection Picks
%         GatherReflectionPicks{ii} = cat(2,ReflectionFBpick{:,ii});
%         Reflection = GatherReflectionPicks{ii};
% %         for jj = 1:length(Radar{ii})
%         parfor (jj = looper, nWorkers)
%             % Cross-Validation for Reflection Velocity Estimation
%             % Multiple Shot Gathers in Population
%             isManyShots = 0;
%             if isManyShots
%                 range = sqrt(([1:length(Radar{ii})]-jj).^2); % Compute Distance
%                 getIx = find(range<=5); % Find Nearby Picks
%                 RefPick = GatherReflectionPicks{ii}(getIx,:) - vertcat(deltaTair{ii,getIx}); % Residual
%                 RefPick = RefPick - vertcat(toAir{ii,getIx})*ones(1,nChan);
%                 stdRefPick = std(RefPick,0,1);
%             end
%             % Single Shot Gather in Population
%             isSingleShot = 1;
%             if isSingleShot
%                 RefPick = Reflection(jj,:) - deltaTair{ii,jj}; % Residual
%                 RefPick = RefPick - toAir{ii,jj};
%             end
% 
%             for kk = 1:100 % 100 Random Draws
%                 % Cross-Validation for Surface Velocity Estimation 7-28-17
%                 nCut = randsample([0,1,2],1);
%                 cutChan = randsample(liveChan,nCut);
%                 xvalChan = liveChan;
%                 cutIx = find(ismember(liveChan,cutChan));
%                 xvalChan(cutIx) = [];                
%                 xvalIx = find(ismember(xvalChan,liveChan));
%                 
%                 xvalRefPick = RefPick(xvalIx);
%                 xvalOffset = offsetArray(xvalIx);
%                 
%                 % Compute Reflected Arrival Velocity and Intercept Time
%                 % OLS Scheme
% %                 [xVref{ii,jj}(kk,1), xToRef{ii,jj}(kk,1), xDepth{ii,jj}(kk,1)] ...
% %                     = Vrms(xvalOffset,xvalRefPick);
%                 % IRLS Scheme
%                 [xVref{ii,jj}(kk,1), xToRef{ii,jj}(kk,1), xDepth{ii,jj}(kk,1)] ...
%                     = VrmsIrls(xvalOffset,xvalRefPick);
%             end
%             
%             % Create Random Sample Population of Surface Density
%             xRhoRef{ii,jj} = DryCrim([xVref{ii,jj}]);
%             
%             % Error Analysis
%             realIx = find(real([xToRef{ii,jj}]));
%             VoRef{ii,jj} = mean([xVref{ii,jj}(realIx)]);
%             VoVarRef{ii,jj} = var([xVref{ii,jj}(realIx)]);
%             toRef{ii,jj} = mean([xToRef{ii,jj}(realIx)]);
%             toVarRef{ii,jj} = var([xToRef{ii,jj}(realIx)]);
%             Depth{ii,jj} = mean([xDepth{ii,jj}(realIx)]);
%             RhoRef{ii,jj} = mean([xRhoRef{ii,jj}(realIx)]);
%             CovZP = cov([xDepth{ii,jj}(realIx)],[xRhoRef{ii,jj}(realIx)]);
%             DepthVar{ii,jj} = CovZP(1,1);
%             RhoRefVar{ii,jj} = CovZP(2,2);
%             CovDepthRho{ii,jj} = CovZP(1,2);
%             
%         end
%         % Concatenate Reflection Travel Time
%         ReflectionTo{ii} = [toRef{ii,:}];
%         ReflectionToVar{ii} = [toVarRef{ii,:}];
%         % Concatenate Reflection Velocity
%         ReflectionVelocity{ii} = [VoRef{ii,:}];
%         ReflectionVelocityVar{ii} = [VoVarRef{ii,:}];
%         % Concatenate Reflector Depth
%         ReflectionDepth{ii} = [Depth{ii,:}];
%         ReflectionDepthVar{ii} = [DepthVar{ii,:}];
%         % Concatenate Reflection Density & Convert Units
%         ReflectionDensity{ii} = [RhoRef{ii,:}];
%         ReflectionDensityVar{ii} = [RhoRefVar{ii,:}];
%         % Concatenate Covariance
%         CovarianceDepthDensity{ii} = [CovDepthRho{ii,:}];
%         
%         % Estimte SWE
%         SWE{ii} = ReflectionDepth{ii}.*ReflectionDensity{ii};
%         % Error Propagation Equation
%         SWEvar{ii} = ReflectionDepthVar{ii}.*ReflectionDensity{ii}.^2 ...
%             + ReflectionDensityVar{ii}.*ReflectionDepth{ii}.^2 ...
%             + 2.*ReflectionDepth{ii}.*ReflectionDensity{ii}.*CovarianceDepthDensity{ii};
%         
%         
%     end
% %     % Smooth Esitmates
% %     Depth = ReflectionDepth{1};
% %     Density = ReflectionDensity{1};
% %     DensityVar = ReflectionDensityVar{1};
% %     DepthVar = ReflectionDepthVar{1};
% %     swe = SWE{1};
% %     sweVar = SWEvar{1};
% %     smoothSWE = conv(swe,hamming(251),'valid')./sum(hamming(251));
% %     smoothSWEvar = conv(sweVar,hamming(251),'valid')./sum(hamming(251));
% %     smoothDensity = conv(Density,hamming(251),'valid')./sum(hamming(251));
% %     smoothDensityVar = conv(DensityVar,hamming(251),'valid')./sum(hamming(251));
% %     smoothDepth = conv(Depth,hamming(251),'valid')./sum(hamming(251));
% %     smoothDepthVar = conv(DepthVar,hamming(251),'valid')./sum(hamming(251));
% %     
% %     % Extrapolate Edges of Convolution
% %     smoothSWE = [ones(1,floor(251/2)).*smoothSWE(1),...
% %     smoothSWE,ones(1,floor(251/2)).*smoothSWE(end)];
% %     smoothSWEvar = [ones(1,floor(251/2)).*smoothSWEvar(1),...
% %     smoothSWEvar,ones(1,floor(251/2)).*smoothSWEvar(end)];
% %     smoothDensity = [ones(1,floor(251/2)).*smoothDensity(1),...
% %     smoothDensity,ones(1,floor(251/2)).*smoothDensity(end)];
% %     smoothDensityVar = [ones(1,floor(251/2)).*smoothDensityVar(1),...
% %     smoothDensityVar,ones(1,floor(251/2)).*smoothDensityVar(end)];
% %     smoothDepth = [ones(1,floor(251/2)).*smoothDepth(1),...
% %     smoothDepth,ones(1,floor(251/2)).*smoothDepth(end)];
% %     smoothDepthVar = [ones(1,floor(251/2)).*smoothDepthVar(1),...
% %     smoothDepthVar,ones(1,floor(251/2)).*smoothDepthVar(end)];
% 
%     % Smooth Esitmates
%     TWT = ReflectionTo{1};
%     TWTvar = ReflectionToVar{1};
%     Depth = ReflectionDepth{1};
%     Density = ReflectionDensity{1};
%     DensityVar = ReflectionDensityVar{1};
%     DepthVar = ReflectionDepthVar{1};
%     swe = SWE{1};
%     sweVar = SWEvar{1};
%     smoothTWT = nonParametricSmooth( 1:length(TWT),TWT,1:length(TWT),251);
%     smoothTWTvar = nonParametricSmooth( 1:length(TWTvar),TWTvar,1:length(TWTvar),251);
%     smoothSWE = nonParametricSmooth( 1:length(swe),swe,1:length(swe),251);
%     smoothSWEvar = nonParametricSmooth( 1:length(sweVar),sweVar,1:length(sweVar),251);
%     smoothDensity = nonParametricSmooth( 1:length(Density),Density,1:length(Density),251);
%     smoothDensityVar = nonParametricSmooth( 1:length(DensityVar),DensityVar,1:length(DensityVar),251);
%     smoothDepth = nonParametricSmooth( 1:length(Depth),Depth,1:length(Depth),251);
%     smoothDepthVar = nonParametricSmooth( 1:length(DepthVar),DepthVar,1:length(DepthVar),251);
%     
%     fprintf('Surface Velocity Analysis Done \n')
%     toc
%     display(' ')
% 
% else    % If No Time Correction is Used ~ Not Recommended
%     for ii = 1:nFiles
%         % Record Time
%         time{ii} = [0:dt:(size(Radar{ii,1},1)-(1+padding(ii)))*dt];
%     end
% end

% %% Create Figure
% if isSWEDISH
% figure();
% subplot(3,1,1)
% shadedErrorBarT8([],smoothSWE,sqrt(smoothSWEvar),1,{'Color',[0.5,0,0],'linewidth',1.5})
% freezeColors
% axis tight
% grid on
% title('Accumulation [m w.e.]') 
% subplot(3,1,2)
% shadedErrorBarT8([],smoothDensity.*1000,sqrt(smoothDensityVar).*1000,1,{'Color',[0,0,.5],'linewidth',1.5})
% freezeColors
% axis tight
% grid on
% title('Density [kgm^{-3}]') 
% subplot(3,1,3)
% shadedErrorBarT8([],smoothDepth,sqrt(smoothDepthVar),1,{'Color',[0,0,0],'linewidth',1.5})
% freezeColors
% axis tight
% grid on
% title('Depth [m]')
% set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
% end
%% Data Export
%     % Create CMP Structure for Export of Data and Coherence Matrix
%     field1 = 'CMP';         value1 = num2cell(selection);
%     field2 = 'Offset';      value2 = selectOffsetArray;
%     field3 = 'Location';    value3 = selectMidPointLoc;
%     field4 = 'Trace';       value4 = selectTraceIxBin;
%     field5 = 'Data';        value5 = selectProCMP;
%     field6 = 'To';          value6 = to;
%     field7 = 'Vo';          value7 = Vo;
%     field8 = 'Time';        value8 = T;
%     field9 = 'Velocity';    value9 = V;
%     field10 = 'Semblance';   value10 = VelX;
%     
%     semblanceCMP = struct(field1,value1,field2,value2,field3,value3,field4,value4, ...
%         field5,value5,field6,value6,field7,value7,field8,value8, field9, value9, field10, value10);
%     
%     % Smash Coherence Matrix
%     stackS = cat(3,VelX{:});
%     stackSemblance = sum(stackS,3);
%     
%     % Normalize Smashed Coherence Matrix
%     SuperSemblance = stackSemblance/max(abs(stackSemblance(:)));
%     % Create CMP Structure for Export of Data and Coherence Matrix
%     field1 = 'DataTrueGain';        value1 = profileCMP{1}{157};
%     field2 = 'DataAGCGain';         value2 = agcCMP{1}{157};
%     field3 = 'Offset';          value3 = offsetArray;
%     field4 = 'RecordTime';        value4 = time{1};
%     field5 = 'StackingTime';      value5 = T;
%     field6 = 'StackingVelocity';    value6 = V;
%     field7 = 'VelocityCoherence';   value7 = veloGram{1}{157};
%     
%     dylanCMP = struct(field1,value1,field2,value2,field3,value3,field4,value4, ...
%         field5,value5,field6,value6,field7,value7);
%     
%     if isWrite
%         cd '/home/tatemeehan/GreenTracs2016/GREENLAND/#ProcessGreenland/GPRprocess/Save'
%         save('TimeCorrectedCMP.mat','dylanCMP','-v7.3')
%         % save('Core7Semblance.txt','V','T','S','-ascii')
%         cd '/home/tatemeehan/GreenTracs2016/GREENLAND/#ProcessGreenland/GPRprocess/TM'
%     end

%     

% 
% % Write Data to .mat File
% if isWrite
%     cd '/home/tatemeehan/GreenTracs2016/GREENLAND/#ProcessGreenland/GPRprocess/Save'
%     save('Core7SemblanceGuideR0CCkm.mat','semblanceCMP','-v7.3')
%     % save('Core7Semblance.txt','V','T','S','-ascii')
%     cd '/home/tatemeehan/GreenTracs2016/GREENLAND/#ProcessGreenland/GPRprocess/TM'
% end
