clear; close all; clc;
%% MxRadar Reads, Sorts, and Processes the Multiplexed Sensors and Software
%  MultiChannel GPR Record. Processed data are stored in the common-offset 
%  domain. Offsets are sorted into common shot gathers for
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
% dataDir = '/home/tatemeehan/GreenTracs2016/GPR_DATA/PulseEKKO/500MHz/6-2-16-Core7-Spur-W';

% GreenTrACS 2017
% dataDir = '/home/tatemeehan/GreenTracs2017/6-19-17-SummitRouteII';
% dataDir = 'E:\GreenTrACS\2017\6-19-17-SummitRouteII';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\6-19-17-SummitRouteII';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\6-18-17-SummitRouteI';
% dataDir = 'E:\GreenTrACS2017\PulseEKKO\500MHz\6-15-17-Core15Core16Traverse';
dataDir = '/home/tatemeehan/GreenTracs2017/6-13-17-Core15SpurW';
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
% dataDir = '/home/tatemeehan/GreenTracs2017/5-31-17-Core11SpurW';
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
% dataDir = '/home/tatemeehan/GreenTracs2017/GPR_DATA/PulseEKKO/1GHz/6-13-17-Core15SpurW';

% Add additional useful pathways
% addpath 'E:\GreenTrACS\2017\GPR_Processing\MultiOffset\Save';
% addpath 'E:\GreenTrACS\2017\GPR_Processing\MultiOffset\SeismicLab\SeismicLab\codes\fx';
% addpath 'D:\GreenTrACS\2017\GPR_Processing\MultiOffset\Save';
% addpath '/run/media/tatemeehan/RED/GreenTrACS/2017/GPR_Processing/MultiOffset/Save'
addpath '/home/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/Save'
addpath '/home/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/Save/matData'
addpath '/home/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/TM'
% addpath '/home/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/SeismicLab/SeismicLab/codes/fx'
% addpath '/sonichome/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/Save'
% addpath '/sonichome/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/SeismicLab/SeismicLab/codes/fx'
addpath '/home/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/time2depth'

%% Control Floe
% Parallel Computing Enabled
isParallel = 1;

% Read Data
isReadSensorsSoftware = 1;     % Read Multiplexed Data
isLoadTimeHorizons = 1;        % Load Previously Picked Time Horizons
isPolarPicker = 0;             % Pick Travel-Time Horizons
isLoadHVA = 0;                 % Load Previous Horizon Velocity Analysis
isLoadMxHL = 0;                % Load Previous MxHL Model Results

% Export Data
isWriteTimeHorizons = 0;

% Process Data
isTrimTWT = 0;          % Truncate Recorded Data for Near-Surface Analysis
isKill = 1;             % Kill Unanted Channels
isMedianSubtraction = 1;% Background Median Subtraction Filter
isFXdecon = 1;          % Fx-Predicitive Deconvolution
isSWEDISH = 1;          % Perform Surface Velocity Analysis
isDepthSection = 1;     % NMO Correction, Stacking, Depth Conversion, Image

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
% load('CalibrationChannelShiftTrough500MHz2017.mat');
load('CalibrationChannelShiftPeak500MHz2017.mat');
% load('CalibrationChannelShiftTrough1GHz2017.mat');
chanShift = Calibration(1).chanShift; % Import chanShift
% Import Cool Colormap
load(['LateNite.mat']);load(['Smoke.mat']); load(['yet_white.mat']);
load(['SplitJet.mat']);

TraverseDistance = [15,15,15];  % Approx. Distance of Radar Files [km]
lineNo = [0,1,2];               % Array of data "LINE" numbers
nFiles = length(lineNo);        % Number of Files
nChan = 9;                      % Number of Recorded Channels
chan =  1:nChan;                % Linear Array of Record Channels
liveChan = chan;

% Establish Tx Rx Geometry for CMP Gathering
nTx = 3;               % Number of Transmitters in Sequence
nRx = 3;               % Number of Receivers in Sequence

% Allocate Memory
time = cell(1,nFiles);
<<<<<<< HEAD

%% Load HVA Results
if isLoadHVA
    load('6-12-17-Core15-Spur-W-HVA-SurfaceForcing-Corrected.mat');
    load('6-12-17-Core15-Spur-W-HVAsmooth-SurfaceForcing-Corrected.mat');
    load('6-12-17-Core15-Spur-W-Bootstraps-Corrected.mat');
end
=======
>>>>>>> 85c68b885e3a89345e5b7f545541f54e66477096
%% Import GPS Information
% Yet Completed
%% Read GPR Data
if isReadSensorsSoftware
<<<<<<< HEAD

    ReadPreProcessRoutine
    
end
%% Concatenate Files
    fprintf('Begin File Concatenation \n')
    tic
    
    catGPR
    
    fprintf('Concatenation Done \n')
    toc
    display(' ')
%% Median Background Subtraction Filter for Reflection Analysis
if isMedianSubtraction
    fprintf('Begin Coherent Noise Removal \n')
    tic
    
    CoherentNoiseRemoval
    
    fprintf('Coherent Noise Removal Done \n')
    toc
    display(' ')
=======
    Rad = cell(1,nFiles);
    trhd = cell(1,nFiles);
    Year = zeros(nFiles,1);
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
        [Rad{ii},hdr1,trhd{ii},dt,f0,~,dx,date] = readSensorsSoftwareData( filepath );
        Year(ii) = str2num(date(1:4));
        % Nominal Frequency GHz
        f0GHz = f0/1000;
        % No. Traces in Multiplexed Data
        [~, multiplexNtrcs] = size(Rad{ii});
        
        % Install Transmitter and Receiver Sequencing and Geometry
        % 500 MHz Offsets
        if f0 == 500
            if nChan == 9
                % 9 Channel Offsets
                txGeo = [0, 0, 0, 1.33, 1.33, 1.33, 2.67, 2.67, 2.67]; % Tx Sequence & Absolute Position
                rxGeo = [4, 8, 12, 4, 8, 12, 4, 8, 12]; % Rx Sequence & Absolute Position
            elseif nChan == 6 && ~isKill
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
            if ii == 1
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
            % Determine Far Offset After Trace Kill
            farOffset = max(offsetArray);           
            farChan = find(offsetArray == farOffset);
            % Determine Near Offset After Trace Kill
            nearOffset = min(offsetArray);          
            nearChan = find(offsetArray == nearOffset);
        end
        
        % Pad Data with Instrument Zero
        if strcmp(date(1:4),'2016')
            padding = 50;
        end
        if strcmp(date(1:4),'2017')
            padding = 0;
        end
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
            reSample = 575 + padding;   % Number of Wanted Samples
            Rad{ii} = Rad{ii}(1:reSample,:);
        end
        
        % Allocation Here
        if ii == 1
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
                
%                 % Sort Trace Indicies
%                 trcIx = traceIx{jj,ii}(:,1:gatherLength - xTrc);
%                 traceIxCell{ii,jj} = trcIx(:);
%                 
%                 % Sort Transmit Position
%                 TxArray{jj,ii} = Array{jj,ii}(:,(1:gatherLength - xTrc)) - txGeo(liveChan(jj));
%                 
%                 % Create SEG-Y Source Location Array
%                 srcxyz{jj,ii} = TxArray{jj,ii};
%                 
%                 % Sort Receive Position
%                 RxMod = length(xArray) - mod(length(xArray),nChan);
%                 RxArray{jj,ii} = xArray(jj:nChan:RxMod) - rxGeo(liveChan(jj));
%                 
%                 % Create SEG-Y Receive Location Array
%                 recxyz{jj,ii} = RxArray{jj,ii};
%                 
%                 % Calculate MidPoint Position (Zero Datum is Tx1 Position)
%                 MidPointX = ((TxArray{jj,ii} + RxArray{jj,ii})/2);
%                 MidPointCell{ii,jj} = MidPointX(:);                 

        end
        fprintf('Signal Processing Done \n')
        toc
        display(' ')
    end
end
%% Concatenate Files
    % Grab Early Time Data For AirWave Cross Correlation Alignment
    % 2016 Approximate Direct Wave Arrival Times
    if strcmp(date(1:4),'2016')
        xcorrWindow = [30, 105, 180, 1, 80, 145, 1, 45, 130;...
                    60, 150, 210, 40, 110, 175, 5, 80, 160] + padding;        
%         xcorrWindow = [30, 105, 180, 0, 80, 145, -20, 45, 130;...
%                     60, 150, 210, 40, 110, 175, 5, 80, 160] + padding;
    end
    % 2017 Approximate Direct Wave Arrival Times
    if strcmp(date(1:4),'2017')
%         xcorrWindow = [175, 250, 85, 140, 230, 75, 100, 180, 25;...
%                     200, 300, 150, 160, 260, 115, 150, 220, 75] + padding;
        xcorrWindow = [90, 155, 220, 65, 125, 195, 45, 110, 175;...
                    115, 190, 260, 90, 165, 265, 75, 135, 205] + padding;
    end
    if isKill
        xcorrWindow(:,killChan) = [];
    end
    
    % Cross-Correlation for multipleFile Alignment Prior to Concatenation
    Radar = xcorrAlignTraces( Radar, xcorrWindow );
    
    % Concatenate Files
    isCat = 1;
    if isCat
        nFiles = 1;
        TraverseDistance = sum(TraverseDistance);
        catRadar = cell(nChan,nFiles);
        for jj = chan            
            catRadar{jj} = [Radar{jj,:}];
        end
         Radar = catRadar; clear catRadar;
    end
    

%% Median Background Subtraction Filter for Reflection Analysis
if isMedianSubtraction
    medR = 1000; % Rank of Median Subtraction Window: W = 2R+1
    % Subtract Median at Later Times
    if f0 == 500
        MedT = xcorrWindow(2,:)+350;
        Taper = 0.25;
        if isTrimTWT
            for ii = 1:nFiles
                for jj = 1:nChan
                    if MedT(jj) > size(Radar{jj,ii},1)
                        MedT(jj) = size(Radar{jj,ii},1);
                    end
                end
            end
        end
%         MedT = 350 + padding;
    end
    
    if f0 == 1000
        MedT = ones(1,nChan);
        Taper = 0.25;
    end
    
    for ii = 1:nFiles
        parfor (jj = chan, nWorkers)
% 
%         for jj = chan;
%             keyboard
            % N samples of topTaper
            topTaper = round((MedT(jj)./2).*Taper);
            % Allocate Unfiltered Signal Matrix
            topRadar = zeros(size(Radar{jj,ii}));
            % Allocate Median Removal Matrix
            medRadar = zeros(size(Radar{jj,ii}));            
            % Solve Percentage for Median Window Tukey Taper
            medTaper = 2.*(topTaper./(size(medRadar,1)-MedT(jj)+topTaper));

            
            
             topRadar(1:MedT(jj),:) = Radar{jj,ii}(1:MedT(jj),:);
             topRadar(1:MedT(jj),:) = topRadar(1:MedT(jj),:).*(tukeywin(MedT(jj),Taper)*ones(1,size(topRadar,2)));
             
             
             medRadar(MedT(jj)-(topTaper-1):end,:) = movingMedianSubtraction...
                 (Radar{jj,ii}(MedT(jj)-(topTaper-1):end,:),medR);
             
             
             tmedL = size(medRadar,1)-MedT(jj)+topTaper;
             medRadar(MedT(jj)-(topTaper-1):end,:).*(tukeywin(tmedL,medTaper)*ones(1,size(medRadar,2)));
%              medRadar(MedT(jj)-round((MedT(jj)./2).*Taper):end,:) = tukeywin(medRadar,Taper);
             Radar{jj,ii} = topRadar+medRadar;
%              Radar{jj,ii}(MedT:end,:) = movingMedianSubtraction(Radar{jj,ii}(MedT:end,:),medR);
%             Radar{jj,ii} = movingMedianSubtraction(Radar{jj,ii},medR);
%             pickRadar{jj,ii} = pickRadar{jj,ii}-median(pickRadar{jj,ii},2)...
%                 *ones(1,size(pickRadar{jj,ii},2));
        end
    end
>>>>>>> 85c68b885e3a89345e5b7f545541f54e66477096
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
%% Branch Processing for Direct and Reflected Arrivals
directRadar = Radar;
%% Snow Water Equivalent and Density for Ice Sheet Height

    if isLoadTimeHorizons
<<<<<<< HEAD
        TimeHorizonFilename = '6-12-17-Core15-Spur-W1-Surface-8chan-TimeHorizon.mat';
        load(TimeHorizonFilename);
        
=======
%         TimeHorizonFilename = '6-2-16-Core7-Spur-W-TimeHorizon.mat';
        TimeHorizonFilename = '6-12-17-Core15-Spur-W1-Surface-8chan-TimeHorizon.mat';
        load(TimeHorizonFilename);
%         if isKill
%         % Remove Killed Offsets
%         DirectFBpick(killChan,:) = [];
%         ReflectionFBpick(killChan,:) = [];
%         end
        % Hacks 10-23-17
        % Remove Noisy Horizon
%         ReflectionFBpick(:,[10]) = [];
%         ReflectionFBpick(:,[3,9,10]) = [];
        % Correct Bad Horizon Pick
%         StaticCorrectionChan = find(offsetArray == 8);
%         ReflectionFBpick{StaticCorrectionChan,6} = ReflectionFBpick{StaticCorrectionChan,6}+4.6;
>>>>>>> 85c68b885e3a89345e5b7f545541f54e66477096
        % Determine Number of Direct Wave Horizons
        nDirectHorizon = size(DirectFBpick,2);
        % Determine Number of Reflection Horizons
        nReflectionHorizon = size(ReflectionFBpick,2);
   
        % Plot Common-offset gathers
        isPlotOffsetGathers = 0;
        isTopMute = 0;
        
<<<<<<< HEAD
        VisualizeDataGathers
        
    end
    
%% Semi-Automatic Radar Wave Picking
=======
        
        if isPlotOffsetGathers
            % AGC Gain for Plotter
            for ii = 1:nFiles
                for jj = chan
                    plotRadar{jj,ii} = AGCgain(Radar{jj,ii},50,2);
                end
            end
%             if isPlotOffsetGathers

            if isTopMute
            % Top Mute for Plotter
            for ii = 1:nFiles
                for jj = chan
                    plotRadar{jj,ii}(1:xcorrWindow(2,jj),:) = ones(xcorrWindow(2,jj),1)*Radar{jj,ii}(1,:);
                end
            end
            end
            
            % Plot Travel-Time Picks
%             time = linspace(-7.45,62.55,350);
%             distance = linspace(0,45,18974);
%             for ii = 1:nFiles
%                 for jj = chan
%                     figure(chan(jj)+(ii-1).*length(chan));...clf;...
%                         imagesc(distance,time,plotRadar{jj,ii}(1:350,:));colormap(Smoke);
%                     title(['Offset: ',num2str(offsetArray(jj)),'m'])
%                     ylabel('Two-Way Time [ns]')
%                     xlabel('Distance [km]')
%                     set(gca,'fontsize',14,'fontweight','bold')
%                     for kk = 1:nDirectHorizon
%                         hold on; plot(distance,...
%                             [DirectFBpick{jj,kk}]-7.45,'k','linewidth',5)
%                     end
%                     for kk = 1:nReflectionHorizon
%                     hold on; plot(distance,...
%                         [ReflectionFBpick{jj,kk}]-7.45,'k','linewidth',5)
%                     end
%                 end
%             end            
            for ii = 1:nFiles
                for jj = chan
                    figure(chan(jj)+(ii-1).*length(chan));clf;...
                        imagesc(plotRadar{jj,ii});colormap(Smoke);
                    title(['Offset: ',num2str(offsetArray(jj)),'m'])
                    for kk = 1:nDirectHorizon
                        hold on; plot(1:length(DirectFBpick{jj,1}),...
                            [DirectFBpick{jj,kk}]./.2,'k','linewidth',2)
                    end
                    for kk = 1:nReflectionHorizon
                    hold on; plot(1:length(ReflectionFBpick{jj,1}),...
                        [ReflectionFBpick{jj,kk}]./.2,'k','linewidth',2)
                    end
                end
            end
            
        end
    end
    
    %% Load HVA Results
    if isLoadHVA
        load('6-12-17-Core15-Spur-W-HVA-SurfaceForcing.mat');
        load('6-12-17-Core15-Spur-W-HVAsmooth-SurfaceForcing.mat');
        load('6-12-17-Core15-Spur-W-Bootstraps.mat');
    end
    %% Semi-Automatic Radar Wave Picking
>>>>>>> 85c68b885e3a89345e5b7f545541f54e66477096
if isSWEDISH
    if isPolarPicker
    fprintf('Begin PolarPicker \n')
    display(' ')
    tic

<<<<<<< HEAD
    PickTravelTimeHorizons    
=======
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
            reflectRadar{jj,ii} = AGCgain(reflectRadar{jj,ii},350,2);
        end
    end
    
    % Plot Common-offset gathers
    isPlotOffsetGathers = 1;
    if isPlotOffsetGathers
        for ii = 1:nFiles
            for jj = chan
                figure(chan(jj)+(ii-1).*length(chan)+100);...
                    imagesc(reflectRadar{jj,ii});colormap(LateNite);
%                 figure(chan(jj)+(ii-1).*length(chan)+200);...
%                     imagesc(directRadar{jj,ii});colormap(LateNite);
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
%     [~, ReflectionFBpick] = polarPickerT8(reflectRadar);
[~, ReflectionFBpick] = polarPickerGetFrame(reflectRadar);
    
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
    
>>>>>>> 85c68b885e3a89345e5b7f545541f54e66477096

    fprintf('PolarPicker Done \n')
    toc
    display(' ')
    end
%% Perform Horizon Velocity Analysis for Estimates of Density, Depth, & SWE
    fprintf('Begin Horizon Velocity Analysis \n')
    display(' ')
    tic
    
    % Toggle Inversion Scheme
    % Air-Coupled Wave Soltion
    isL1Air = 0;
    isL2Air = 1;
    % Surface-Coupled Wave
    isL1LMO = 0;
    isL2LMO = 1;
    % Reflected Wave
    isL1NMO = 0;
    isL2NMO = 1;
    
    % Run Memory Allocation
    HVAmemoryAllocation
    
    % Bias Calibration at Core 15W Pit (Observed - Estimated);
    % This is a 1.72% Adjustment
    velocityBias = 0.004;      % [m/ns]
    densityBias = -.0258902;    % [g/cm3]
    depthBias = 0.7647;          % [m] 2016 Horizon
    
<<<<<<< HEAD
    HorizonVelocityAnalysis


end   
%% Apply Residual Trace Shifts to Data Gathers
if isSWEDISH || isLoadHVA && ~isLoadMxHL

    DataStaticShift
=======


% Extract AirWave Picks and Perform Residual Subtraction Velocity Analysis
%     parfor (ii = 1:nFiles, nWorkers - (nWorkers-nFiles)) 
    for ii = 1:nFiles
%         looper = 1:size(Radar{ii},2);
          looper = 1:nTrace;
        % Concatenate Air Wave Picks
        for dh = 1:nDirectHorizon
            % Concatenate Primary Reflection Picks for Horizon hh
            GatherDirectPicks{dh,ii} = cat(2,DirectFBpick{:,dh,ii});
            DirectPicks = GatherDirectPicks{dh,ii};
            % Air Wave Velocity Analysis
            if dh  == 1
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
                        % IRLS Scheme
                        if isL1Air
                            [xVdir{dh,jj,ii}(kk,1), xToDir{dh,jj,ii}(kk,1)] = DirectWaveIrls(xvalOffset,xvalAirPick);
                        end
                        % OLS Scheme
                        if isL2Air
                            [xVdir{dh,jj,ii}(kk,1), xToDir{dh,jj,ii}(kk,1)] = DirectWave(xvalOffset,xvalAirPick);
                        end
                       
                    end
                    
                    % Residual Error Analysis
                    VoAir = mean([xVdir{dh,jj,ii}]);
                    tmpTo = mean([xToDir{dh,jj,ii}]);%keyboard
%                     
%                     residualOffset = offsetArray - ((AirPick-tmpTo).*VoAir);
%                     residualTime = residualOffset./VoAir;
%                     residualTime = residualOffset./.299;
                    
                    % Find Velocity Residual for Additional Trace Shifts & Calulate Residual Time
%                     keyboard
%                     deltaT{ii,jj} = (offsetArray./VoAir) - (offsetArray./0.299);
%                     ResAirPick = AirPick - deltaT{ii,jj};
% ResAirPick = AirPick - residualTime;
                    % Fix Residual Corrections 11-3-17
                    deltaT{ii,jj} = (AirPick-tmpTo) - (offsetArray./0.299);
                    ResAirPick = AirPick - deltaT{ii,jj};
%                     keyboard
                    % Recalculate Airwave Moveout Velocity
                    [VoDir{dh,jj,ii}, toDir{dh,jj,ii}] = DirectWaveIrls(offsetArray,ResAirPick);
%                     [VoDir{dh,jj,ii}, toDir{dh,jj,ii}] = DirectWave(offsetArray,ResAirPick);

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
                        ranger = sqrt(([1:nTrace]-jj).^2); % Compute Distance
%                         ranger = sqrt(([1:size(Radar{ii},2)]-jj).^2); % Compute Distance
                        getIx = find(ranger<=shotRange); % Find Nearby Picks
                        DirPick = DirectPicks(getIx,:) - vertcat(deltaT{ii,getIx}); % Residual Static Corection
%                         DirPick = DirPick - vertcat(AirTo{ii,getIx})*ones(1,nChan); % Time-Zero Static Shift
                        pickPool = size(DirPick,1);
                        xvalPool = 1:pickPool;
                        % Cross-Validation for Surface Velocity Estimation 10-12-17
                        for kk = 1:1250 % 1250 Random Draws
                            nCut = randsample([0,1],1);
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
                            if isL2LMO
                            % OLS Scheme
                            [xVdir{dh,jj,ii}(kk,1), xToDir{dh,jj,ii}(kk,1)]...
                              = DirectWave(xvalOffset,xvalDirPick);
                            end
                            % IRLS Scheme
                            if isL1LMO
                            [xVdir{dh,jj,ii}(kk,1), xToDir{dh,jj,ii}(kk,1)]...
                                = DirectWaveIrls(xvalOffset,xvalDirPick);
                            end
                        end
                    end
                    % Single Shot Gather in Population
                    isSingleShot = 0;
                    if isSingleShot
                        % Re-Extract Direct Wave Picks after Residual Shifts
                        DirPick = DirectPicks(jj,:) - deltaT{ii,jj}; % Residual Static Corection
%                         DirPick = DirPick - AirTo{ii,jj}*ones(1,nChan); % Time-Zero Static Shift                       
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
                            if isL2LMO
                            % OLS Scheme
                            [xVdir{dh,jj,ii}(kk,1), xToDir{dh,jj,ii}(kk,1)]...
                              = DirectWave(xvalOffset,xvalDirPick);
                            end
                            % IRLS Scheme
                            if isL1LMO
                            [xVdir{dh,jj,ii}(kk,1), xToDir{dh,jj,ii}(kk,1)]...
                                = DirectWaveIrls(xvalOffset,xvalDirPick);
                            end
                        end
                    end
                    % Correct Residual Direct Arrival Intercept Time
                    if dh == 2
                        ResidualToDir = AirTo{ii,jj} - mean(xToDir{dh,jj,ii});
                        AirTo{ii,jj} = AirTo{ii,jj} - ResidualToDir;
                    end
                    % Apply To Static Shift Post Residual Correction
                    xToDir{dh,jj,ii} = xToDir{dh,jj,ii} - AirTo{ii,jj};
                    % Apply Velocity Bias Correction Factor 1.72%
                    xVdir{dh,jj,ii} = xVdir{dh,jj,ii} + velocityBias;
                    % Estimate Wavelet Depth
                    % Sounding Depth at One Half Wavelength
                    waveLength = (xVdir{dh,jj,ii}./f0GHz)./2;
                        % L1 Norm Forces Variance of To Positive Downward
                    xDepthDir{dh,jj,ii} = xVdir{dh,jj,ii}.*abs(xToDir{dh,jj,ii})...
                        + waveLength;
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
                % Apply Residual Corrections to AirWave Arrival Time
                if dh == 2
                    toDir(1,:,ii) = AirTo(ii,:);
                    DirectTo{1,ii} = [toDir{1,:,ii}] - [AirTo{ii,:}];
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
%             if dh > 1
%                 % Apply Core15 Spur W Bias Correction
%                 DirectDensity{dh,ii} = [RhoDir{dh,:,ii}]+densityBias;
%             else
                DirectDensity{dh,ii} = [RhoDir{dh,:,ii}];%+densityBias;
%             end
            
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
        looper = 1:nTrace;
%         looper = 1:size(Radar{ii},2);
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
                ranger = sqrt(([1:nTrace]-jj).^2); % Compute Distance   
%                 keyboard
%                 ranger = sqrt(([1:size(Radar{ii},2)]-jj).^2); % Compute Distance
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
                    nCut = randsample([0,1],1);
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
                    if isL1NMO
                        [xVref{rh,jj,ii}(kk,1), xToRef{rh,jj,ii}(kk,1), xDepth{rh,jj,ii}(kk,1)] ...
                            = VrmsIrls(xvalOffset,xvalRefPick);
                    end
                    % OLS Scheme
                    if isL2NMO
                        [xVref{rh,jj,ii}(kk,1), xToRef{rh,jj,ii}(kk,1), xDepth{rh,jj,ii}(kk,1)] ...
                            = Vrms(xvalOffset,xvalRefPick);
                    end
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
                    % IRLS Scheme
                    if isL1NMO
                        [xVref{rh,jj,ii}(kk,1), xToRef{rh,jj,ii}(kk,1), xDepth{rh,jj,ii}(kk,1)] ...
                            = VrmsIrls(xvalOffset,xvalRefPick);
                    end
                    % OLS Scheme
                    if isL2NMO
                        [xVref{rh,jj,ii}(kk,1), xToRef{rh,jj,ii}(kk,1), xDepth{rh,jj,ii}(kk,1)] ...
                            = Vrms(xvalOffset,xvalRefPick);
                    end
                    end
                end
            % Apply Velocity Bias Correction Factor % 1.72
            xVref{rh,jj,ii} = xVref{rh,jj,ii} + velocityBias;
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
        % Apply Bias Correction
        ReflectionDepth{rh,ii} = [HorizonDepth{rh,:,ii}].*(depthBias);
        ReflectionDepthVar{rh,ii} = [HorizonDepthVar{rh,:,ii}];
        %Concatenate Layer Thicnkess
        IntervalThickness{rh,ii} = [HoInt{rh,:,ii}];
        IntervalThicknessVar{rh,ii} = [HoIntVar{rh,:,ii}];
        % Concatenate Reflection Density
        % % Apply Core 15 Spur W Bias Correction
        ReflectionDensity{rh,ii} = [RhoRef{rh,:,ii}];%+densityBias;
        ReflectionDensityVar{rh,ii} = [RhoRefVar{rh,:,ii}];
        % Concatenate Interval Density
        IntervalDensity{rh,ii} = [RhoInt{rh,:,ii}];
        IntervalDensityVar{rh,ii} = [RhoIntVar{rh,:,ii}];
        % Concatenate Covariance
        CovarianceDepthDensity{rh,ii} = [CovDepthRho{rh,:,ii}];
        CovarianceThicknessDensity{rh,ii} = [CovThicknessRho{rh,:,ii}];
        
        % Estimte SWE
        % Average over 1.5 years Winter 2016 - Summer 2017
        SWE{rh,ii} = ReflectionDepth{rh,ii}.*ReflectionDensity{rh,ii}./1.5;
        SWEint{rh,ii} = IntervalThickness{rh,ii}.*IntervalDensity{rh,ii}./1.5;
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
%         if nDirectHorizon > 1;
%             % Constrained Inversion
%             % Concatenate Reflection Picks
%             groupReflectionTravelTimes = reshape(cat(3,ReflectionFBpick{:,:,ii}),size(ReflectionFBpick{1,1,ii},1),nChan,nReflectionHorizon);
%             groupReflectionTravelTimes = permute(groupReflectionTravelTimes,[2,3,1]);
%             groupDirectTravelTimes = reshape(cat(3,DirectFBpick{:,:,ii}),size(DirectFBpick{1,1,ii},1),nChan,nDirectHorizon);
%             groupDirectTravelTimes = permute(groupDirectTravelTimes,[2,3,1]);
%             parfor (jj = looper, nWorkers)
% %             for jj = looper
%                 [VrmsProfile{ii,jj}, ToProfile{ii,jj}, ZoProfile{ii,jj},...
%                     VintProfile{ii,jj}, HintProfile{ii,jj}]...
%                     = constrainVelocities( groupReflectionTravelTimes(:,:,jj), AirTo{ii,jj}, ...
%                     xToDir(2,jj), xToRef(:,jj), xVdir(2,jj), xVref(:,jj), xDepthDir(2,jj),...
%                     xDepth(:,jj), offsetArray);
%             end
%             Vdix = cell2mat(VintProfile);
%             VRMS = cell2mat(VrmsProfile);
%             Zz = cell2mat(ZoProfile);
%         end

    end
end   
%% Apply Residual Trace Shifts
if isSWEDISH || isLoadHVA && ~isLoadMxHL
% Apply Residual Channel Shifts to Travel Time Picks
% Allocate Memory For Travel-Time Export
exportReflectionTravelTimes = cell(nFiles,1);
exportDirectTravelTimes = cell(nFiles,1);

for ii = 1:nFiles
    % Re-Package Travel-Time Picks
    groupReflectionTravelTimes = reshape(cat(3,ReflectionFBpick{:,:,ii}),size(ReflectionFBpick{1,1,ii},1),nChan,nReflectionHorizon);
    exportReflectionTravelTimes{ii} = permute(groupReflectionTravelTimes,[2,3,1]);
    groupDirectTravelTimes = reshape(cat(3,DirectFBpick{:,:,ii}),size(DirectFBpick{1,1,ii},1),nChan,nDirectHorizon);
    exportDirectTravelTimes{ii} = permute(groupDirectTravelTimes,[2,3,1]);
    for jj = 1:size(deltaT,2)
%         keyboard
        % Apply Residual and Static Corrections for Direct Waves
        exportDirectTravelTimes{ii}(:,:,jj) = exportDirectTravelTimes{ii}(:,:,jj) ...
            - ((ones(nDirectHorizon,1)*deltaT{ii,jj})'...
            + (ones(nDirectHorizon,1)*(ones(1,nChan).*AirTo{ii,jj}))');
        % Apply Residual and Static Corrections for Reflected Waves
        exportReflectionTravelTimes{ii}(:,:,jj) = exportReflectionTravelTimes{ii}(:,:,jj) ...
            - ((ones(nReflectionHorizon,1)*deltaT{ii,jj})'...
            + (ones(nReflectionHorizon,1)*(ones(1,nChan).*AirTo{ii,jj}))');
    end
end

                
% Apply Residual and Static Shifts to Data Amplitudes
% dumRadar = Radar;
for ii = 1:nFiles
    parfor (jj = 1:nChan, nWorkers)
%         dumRadar{jj,ii} = Radar{jj,ii}(1251:end,:);
%     for jj = 1:nChan
        % Remove AirWave 0verhead Time
        Overhead = Radar{jj,ii};
        OverheadIx = round(AirTo{ii,jj}./dt);
        Radar{jj,ii} = [Overhead(OverheadIx+1:end,:);Overhead(1:OverheadIx,:)];
%         dumRadar{jj,ii}(OverheadIx,:) = 1;%[];
% %         Radar{jj,ii}(1:round(AirTo{ii,jj}./dt),:) = [];
%         Radar{jj,ii} = [dumRadar{jj,ii};Overhead];
        % Create Temporary Dimension Handle
        tmpLength = size(Radar{jj,ii},1);
        noise = mean(Radar{jj,ii}(1,:));
        directNoise = mean(directRadar{jj,ii}(1,:));
        for kk = 1:size(deltaT,2)
%             keyboard
            % Apply Systematic Channel Shifts
            shiftSample = round(deltaT{ii,kk}(jj)./dt);
            
            if shiftSample < 0
                % Direct Arrival Data
                tmp = padarray(directRadar{jj,ii}(:,kk),[abs(shiftSample),0],...
                    directNoise,'pre');
                directRadar{jj,ii}(:,kk) = tmp(1:tmpLength,:);
                % Reflection Data
                tmp = padarray(Radar{jj,ii}(:,kk),[abs(shiftSample),0],...
                    noise,'pre');
                Radar{jj,ii}(:,kk) = tmp(1:tmpLength,:);
                
            elseif shiftSample > 0
                % Direct Wave Data
                tmp = padarray(directRadar{jj,ii}(:,kk),[abs(shiftSample),0],...
                    directNoise,'post');
                directRadar{jj,ii}(:,kk) = tmp(abs(shiftSample)+1:end,:);
                % Reflection Data
                tmp = padarray(Radar{jj,ii}(:,kk),[abs(shiftSample),0],...
                    noise,'post');
                Radar{jj,ii}(:,kk) = tmp(abs(shiftSample)+1:end,:);
                
            else
                directRadar{jj,ii}(:,kk) = directRadar{jj,ii}(:,kk);
                Radar{jj,ii}(:,kk) = Radar{jj,ii}(:,kk);
                
            end
        end
    end
end    
>>>>>>> 85c68b885e3a89345e5b7f545541f54e66477096
    
    fprintf('Horizon Velocity Analysis Done \n')
    toc
    display(' ')
<<<<<<< HEAD
=======
% 
% else    % If No Time Correction is Used ~ Not Recommended
%     for ii = 1:nFiles
%         % Record Time
%         time{ii} = [0:dt:(size(Radar{ii,1},1)-(1+padding(ii)))*dt];
%     end
>>>>>>> 85c68b885e3a89345e5b7f545541f54e66477096
end
%% Extrapolate Surface Density from Radar Forcing
% Bias Calibration at Core 15W Pit (Observed - Estimated);
% Applied Within HVA

<<<<<<< HEAD
    MxHLforcingParameter

=======
if ~isSWEDISH
    % Determine Maximum File Size
    nTrace = zeros(1,nFiles);
    for ii = 1:nFiles
        nTrace(ii) = size(Radar{1,ii},2);
    end
    % Maxiumum Traces to Loop Over
    nTrace = max(nTrace);
end

if ~isLoadHVA || ~isLoadMxHL
% Allocation
looper = 1:nTrace;
SurfaceDensity = cell(nFiles,1);
SurfaceDensityVar = cell(nFiles,1);
SurfaceDensification = cell(nFiles,1);
SurfaceDensificationVar = cell(nFiles,1);
ForcingDensity = cell(nFiles,1);
ForcingDensityVar = cell(nFiles,1);
ForcingDepth = cell(nFiles,1);
ForcingDepthVar = cell(nFiles,1);

for ii = 1:nFiles
    MCsamples = 250;
    % Allocation
    surfVar = zeros(length(looper),1);
    rateVar = zeros(length(looper),1);
    surfDensity = zeros(length(looper),1);
    densityRate = zeros(length(looper),1);
    
    % MC BootStrapping Regression for Surface Density Extrapolation
    parfor (jj = looper, nWorkers)
        % Ensure Real Indicies
        realIx = find(real([xToRef{1,jj,ii}]));
        % Create Indicies for Sampling
        dirIx = 1:length([xRhoDir{2,jj,ii}]);
        refIx = 1:length([xRhoRef{1,jj,ii}(realIx)]);
        
        % Monte Carlo Sampling Depth and Density
        dirSample = randsample(dirIx,MCsamples,'true');
        refSample = randsample(refIx,MCsamples,'true');
        
        % dh = 2 Surface Wave Horizon
        xStrapDepthDir = xDepthDir{2,jj,ii}(dirSample);
        xStrapRhoDir = xRhoDir{2,jj,ii}(dirSample);

        % rh = 1 Primary Reflection Horizon
        xStrapDepthRef = xDepth{1,jj,ii}(refSample);
        xStrapRhoRef = xRhoRef{1,jj,ii}(refSample);
        
        % Allocation
        surfRho = zeros(MCsamples,1);
        rhoRate = zeros(MCsamples,1);
        for kk = 1:250
            % L2 Norm Regression
            G = [1,xStrapDepthDir(kk);1,xStrapDepthRef(kk)];
            d = [xStrapRhoDir(kk);xStrapRhoRef(kk)];
            m = G\d;
            surfRho(kk) = m(1);
            rhoRate(kk) = m(2);
        end
        % Estimate BootStrapped Sample Variance
        surfVar(jj) = var(surfRho);
        rateVar(jj) = var(rhoRate);
        surfDensity(jj) = mean(surfRho);
        densityRate(jj) = mean(rhoRate);
    end
    SurfaceDensity{ii} = nonParametricSmooth( 1:length(surfDensity),...
        surfDensity,1:length(surfDensity),smoothR);
    SurfaceDensityVar{ii} = nonParametricSmooth( 1:length(surfVar),...
        surfVar,1:length(surfVar),smoothR);
    SurfaceDensification{ii} = nonParametricSmooth( 1:length(densityRate),...
        densityRate,1:length(densityRate),smoothR);
    SurfaceDensificationVar{ii} = nonParametricSmooth( 1:length(rateVar),...
        rateVar,1:length(rateVar),smoothR);
%     SurfaceDensity{ii} = surfDensity;
%     SurfaceDensityVar{ii} = surfVar;
%     SurfaceDensification{ii} = densityRate;
%     SurfaceDensificationVar{ii} = rateVar;

% The Forcing Density of Herron-Langway (1980) is to Good Approximation the
% Average of the two Radar Surface Estimates.
ForcingDensity{ii} = mean([dhDensity{2,ii},Density{1,ii}],2);
ForcingDensityVar{ii} = mean([dhDensityVar{2,ii},DensityVar{1,ii}],2);
ForcingDepth{ii} = mean([dhDepth{2,ii},Depth{1,ii}],2);
ForcingDepthVar{ii} = mean([dhDepthVar{2,ii},DepthVar{1,ii}],2);
    
end
end
>>>>>>> 85c68b885e3a89345e5b7f545541f54e66477096
%% Create Depth Section and Density Model
if isDepthSection && ~isLoadMxHL
    %% Stacking Velocity Model
    fprintf('Begin Stacking Velocity Extrapolation \n')
    tic
<<<<<<< HEAD
    
    VelocityModel
=======
    % Allocation
    AverageAccumulation = cell(nFiles,1);
    StackingVelocity = cell(nFiles,1);
    StackingDensity  = cell(nFiles,1);
    StackingTime = cell(nFiles,1);
    HerronLangwayDensity = cell(nFiles,1);
    HerronLangwayAge = cell(nFiles,1);
    Traverse = cell(nFiles,1);
    TraverseX = cell(nFiles,1);
    StackDepth = cell(nFiles,1);
    
    % Annual Accumulation Correction from Pit 15W [.Percent]
    % This is the depth Bias Correction Factor
%     normSWE = 0.0990;
    % Annual Accumulation Correction from Core 15 Chemistry [mwe]
    % 50 year mean is ~0.299 [mwe]
    accumulationAvg = 0.299;


    % Estimate of Mean Annual Temperature
%     annualT = -18; % [C]     
%     annualT = -19; % [C]    
%     annualT = -20; % [C]
    annualT = -21.5; % [C]    

    % Array of Model Depths
    maxDepth = 30.001;
    nDepth = 301;
    StackZ = (linspace(0.001,maxDepth,nDepth))';
    % Interval Thickness
    StackH = diff(StackZ);StackH = [StackZ(1);StackH];
    
    for ii = 1:nFiles
        % Array of Approximate Distance
        Traverse{ii} = linspace(0,TraverseDistance(ii),length(dhDensity{2,1}));
        % Distance Grid
        TraverseX{ii} = ones(length(StackZ),1)*Traverse{ii};
        % Depth Grid
        StackDepth{ii} = StackZ*ones(1,length(Traverse{ii}));
        for jj = 1:length(ForcingDensity{ii})
            % Extract Radar Forcing for Herron-Langway Model
            surfHL = ForcingDensity{ii}(jj);
%             accuHL = SnowWaterEqv{1,ii}(jj);% - (normSWE.*SnowWaterEqv{1,ii}(jj));
            
            accuHL = accumulationAvg + (SnowWaterEqv{1,ii}(jj) - mean(SnowWaterEqv{1,ii}(1:100)));
            AverageAccumulation{ii}(jj) = accuHL;
            % Impose Herron-Langway (1980) Density Model
            [HerronLangwayDensity{ii}(:,jj),HerronLangwayAge{ii}(:,jj)] = herronLangway(StackZ,annualT,surfHL,accuHL);
            
            %%% Apply Surface Correction to Herron-Langway Model Here %%%
            radHLdepth0 = StackZ(2);
            radHLrho0 = SurfaceDensity{ii}(jj);
            [~, radHLdepthIx1] = min(abs(StackZ - dhDepth{2,ii}(jj)));
            radHLdepth1 = StackZ(radHLdepthIx1);
            radHLrho1 = dhDensity{2,ii}(jj);
            [~, radHLdepthIx2] = min(abs(StackZ - ForcingDepth{ii}(jj)));
            radHLdepth2 = StackZ(radHLdepthIx2);
            radHLrho2 = ForcingDensity{ii}(jj);
            [~, radHLdepthIx3] = min(abs(StackZ - Depth{1,ii}(jj)));
            radHLdepth3 = StackZ(radHLdepthIx3);
            radHLrho3 = Density{1,ii}(jj);
            
            % L2 Norm Regression
            G = [1,radHLdepth0; 1,radHLdepth1;1,radHLdepth2; 1,radHLdepth3];
            d = [radHLrho0;radHLrho1;radHLrho2;radHLrho3];
            m = G\d;
            surfRhoHL = m(1);
            rhoRateHL = m(2);
            
            %%% Average Radar Derived Surface Density %%%
            avgSurfacePiece = rhoRateHL.*StackZ(1:radHLdepthIx3) + surfRhoHL;
           
            %%% Causal Integration Inversion for Surface Density %%%
            G = tril(ones(length(1:radHLdepthIx3)));
            d = avgSurfacePiece.*sum(G,2);
            m = G\d;
            SurfacePiece = m;
            
            % Find Intersection
            [~, PieceIx] = min(abs(HerronLangwayDensity{ii}(1:length(SurfacePiece),jj)-SurfacePiece));
            % Append Surface Piece
            HerronLangwayDensity{ii}(1:PieceIx-1,jj) = SurfacePiece(1:PieceIx-1);
            
            % Compute Weights for Cumulative Average Density
            stackW = StackZ.^(-1).*tril(ones(length(StackH),1)*StackH');

            % Cumulative Average Density is a Proxy for Stacking Velocity
            StackingDensity{ii}(:,jj) = stackW*HerronLangwayDensity{ii}(:,jj);
            
            % Estimate Stacking Velocity in Time
            StackingVelocity{ii}(:,jj) = DryCrimVRMS(StackingDensity{ii}(:,jj));
            StackingTime{ii}(:,jj) = 2.*StackZ./StackingVelocity{ii}(:,jj);
                      
        end
    end
>>>>>>> 85c68b885e3a89345e5b7f545541f54e66477096
    
    fprintf('Velocity & Density Extrapolation Done \n')
    toc
    display(' ')
    

    %% Normal Moveout Correction
    fprintf('Begin Normal Moveout Correction \n')
    tic
<<<<<<< HEAD
    
    NormalMoveoutCorrection    
    
=======
    Stretch = 100;
    % Allocation
    xStack = cell(nFiles,1);
    tStack = cell(nFiles,1);
    vStack = cell(nFiles,1);
    zStack = cell(nFiles,1);
    DepthAxis = cell(nFiles,1);
    DepthMatrix = cell(nFiles,1);
    RadarNMO = Radar;
    for ii = 1:nFiles
        jj = 1;
        isInterpolate = 1;
        % Perform NMO Correction with Grid Interpolation
        [RadarNMO{jj,ii},xStack{ii},tStack{ii},vStack{ii},~] = ...
            commonOffsetNMO(Radar{jj,ii},dt,f0,offsetArray(jj),TraverseX{ii},StackingTime{ii},StackingVelocity{ii},Stretch,isInterpolate);
        % Grid of Depths
        zStack{ii} = (vStack{ii}.*tStack{ii})./2;
        % Axis for Depth Image
        DepthAxis{ii} = (linspace(min(zStack{ii}(:)),max(zStack{ii}(:)),size(RadarNMO{jj,ii},1)))';
        DepthMatrix{ii} = DepthAxis{ii}(:)*ones(1,size(RadarNMO{jj,ii},2));
        
        % Allocation for Parfor Overhead
        dumX = xStack{ii};
        dumT = tStack{ii};
        dumV = vStack{ii};
        parfor (jj = 2:nChan, nWorkers)
            % Perform NMO Correction on Remaining Channels w/o Interpolation
            isInterpolate = 0;
            [RadarNMO{jj,ii},~,~,~,~] = ...
                commonOffsetNMO(Radar{jj,ii},dt,f0,offsetArray(jj),dumX,dumT,dumV,Stretch,isInterpolate);
        end
        clear('dumX','dumT','dumV');
    end
>>>>>>> 85c68b885e3a89345e5b7f545541f54e66477096
    fprintf('Normal Moveout Correction Done \n')
    toc
    display(' ')
    %% Stack Gathers
    fprintf('Begin Offset-Gather Stacking \n')
    tic
<<<<<<< HEAD
    
    StackGathers
=======
    % Partial Stack
    % Selected Offsets 4 and 12 m for Stack. The Late Times of 4m are muted
    RadarStack = cell(nFiles,1);
    [~,ix4] = min(abs(offsetArray - 4));
    [~,ix12] = min(abs(offsetArray - 12));
    for ii = 1:nFiles
        % Bottom Mute
        stackTaper = zeros(size(RadarNMO{ix4,ii},1),1);
        stackTaper(1:200,1) = tukeywin(size(RadarNMO{ix4,ii}(1:200,:),1),0.05);
        stackTaper = stackTaper*ones(1,size(RadarNMO{ix4,ii},2));
        % Stack Offsets
        RadarStack{ii} = (RadarNMO{ix4,ii}.*stackTaper)+RadarNMO{ix12,ii};
        clear stackTaper
    end
%     clear RadarNMO
>>>>>>> 85c68b885e3a89345e5b7f545541f54e66477096
    
    fprintf('Stacking Done \n')
    toc
    display(' ')
    %% Time to Depth Conversion
    fprintf('Begin Time-Depth Conversion \n')
    tic
<<<<<<< HEAD
    
    TimeToDepthConversion
=======
    % Allocation
    RadarDepth = RadarStack;  
    
    for ii = 1:nFiles
        RadDepth = zeros(size(RadarStack{ii}));
        RadStack = RadarStack{ii};
        DepthAxe = DepthAxis{ii};
        zStak = zStack{ii};
        tStak = tStack{ii};
        parfor (kk = 1:size(RadarStack{ii},2), nWorkers)
            % [Time, Depth] - Array to be resampled
            tzcurve = [tStak(:,kk),zStak(:,kk)];
            % Perform Depth Conversion
            [RadDepth(:,kk)] = timeDepthConversion(RadStack(:,kk),tzcurve,DepthAxe);
        end
        RadarDepth{ii} = RadDepth;
    end
    clear('RadarStack','RadDepth','RadStack','DepthAxe','zStak','tStak');
>>>>>>> 85c68b885e3a89345e5b7f545541f54e66477096
    
    fprintf('Time-Depth Conversion Done \n')
    toc
    display(' ')
    %% Post-Stack FX-Deconvolution
    if isFXdecon
        fprintf('Begin Post-Stack FX-Predictive Deconvolution \n')
        tic
        
        % Select Temporal or Spatial Deconvolution Gates: 1 on; 0 off;
        isTemporalFXdecon = 1;
        isSpatialFXdecon = 0;
        
        PostStackFXdeconTX
        
        fprintf('Post-Stack FX-Deconvolution Done \n')
        toc
        display(' ')
    end
    %% Resize Density Model
    fprintf('Begin Interpolation of MxHL Model \n')
    tic
<<<<<<< HEAD

    MxHLmodel

=======
    % Interpolate Density Model to Resolution of Radar Depth Section
    AgeModel = cell(nFiles,1);
    Ages = cell(nFiles,1);
    DensityModel = cell(nFiles,1);
    DensityAnomalyModel = cell(nFiles,1);
    AvgDensityModel = cell(nFiles,1);
    MeanDensityDeviation = cell(nFiles,1);
    SurfaceDensityDeviation = cell(nFiles,1);
    DepthAge = cell(nFiles,1);
    AvgAccumulation = cell(nFiles,1);
    MxAccumulationRate = cell(nFiles,1);
    MxHLAccumulationRate = cell(nFiles,1);
    ix550 = zeros(1,size(DensityModel{ii},2));
    depth550 = zeros(1,size(DensityModel{ii},2));
    for ii = 1:nFiles
        % Resize Depth-Density Model
        DensityModel{ii} = griddata(TraverseX{ii},StackDepth{ii},...
            HerronLangwayDensity{ii},xStack{ii},DepthMatrix{ii},'linear');
        % Fill in NaN with Nearest Extrapolation
        DensityModel{ii}(1,:) = DensityModel{ii}(2,:);
        
        % Resize Depth-AvgDensity Model
        AvgDensityModel{ii} = griddata(TraverseX{ii},StackDepth{ii},...
            StackingDensity{ii},xStack{ii},DepthMatrix{ii},'linear');
        % Fill in NaN with Nearest Extrapolation
        AvgDensityModel{ii}(1,:) = AvgDensityModel{ii}(2,:);        
        
        % Resize Depth-Age Model
        AgeModel{ii} = griddata(TraverseX{ii},StackDepth{ii},...
            HerronLangwayAge{ii},xStack{ii},DepthMatrix{ii},'linear');
        % Fill in NaN with Zero Age
        AgeModel{ii}(1,:) = 0;
        % Determine Critical Depth 
        for kk = 1:size(DensityModel{ii},2)
            ix550(kk) = find(DensityModel{ii}(:,kk)>.550,1);
            depth550(kk) = DepthMatrix{ii}(ix550(kk),kk);
        end
        
        % Density Anomaly from the Mean
        DensityAnomalyModel{ii} = DensityModel{ii} - mean(DensityModel{ii},2);

        MeanDensityDeviation{ii} = mean(DensityAnomalyModel{ii},1);
        SurfaceDensityDeviation{ii} = DensityAnomalyModel{ii}(1,:);
%         DensityDeviation2D{ii} = DensityAnomalyModel{ii}-ones(size(DensityAnomalyModel,1))*MeanDensityDeviation{ii};

        % Calculate Accumulation Rate from Maximum Isochrone to Surface
        
        MaxAge = floor(min(AgeModel{ii}(end,:)));
%         Ages = [1,3,7,12,MaxAge];
        Ages{ii} = [1,3:3:MaxAge];
%         Isochrones = abs(Year(ii) - [1,3,7,12,MaxAge]);
        Isochrones = abs(Year(ii) - Ages{ii});
        for kk = 1:length(Isochrones)
        [~, AgeIx] = min(abs(AgeModel{ii}-Ages{ii}(kk)));
        AgeIx = sub2ind(size(DepthMatrix{ii}),AgeIx,1:length(AgeIx));
        DepthAge{ii}(:,kk) = DepthMatrix{ii}(AgeIx');
        AvgAccumulation{ii}(:,kk) = AvgDensityModel{ii}(AgeIx').*DepthAge{ii}(:,kk);
        end
        % This Calculation is Redundant Due Circularity of Methodology
        MxHLAccumulationRate{ii} = AvgAccumulation{ii}(:,find(Ages{ii} == MaxAge))./MaxAge;
        % MxRadar Evaluated Contemporary Accumulation Rate
        % This result is Bias Corrected (Depth and Density) is Final
        MxAccumulationRate{ii} = SnowWaterEqv{1,ii};% - (normSWE.*SnowWaterEqv{1,ii});
        % Uncertainties Are Developed During HVA Bootstrapping
        Isochrones = DepthAge;


    end
>>>>>>> 85c68b885e3a89345e5b7f545541f54e66477096
    fprintf('MxHL Model Interpolation Done \n')
    toc
    display(' ')
    end
    %% Plot Wiggle OverLay Velocity
if isLoadMxHL
    MxHLFilename = 'GTC15SpurWMxHL.mat';
    % Load MxHL structure
    load(MxHLFilename);
    % Write Unpacked Structure Variables to .m file
    structvars(GTC15SpurWMxHL);
    % Run .m to Evaluate structvars
    tmpvars
    delete tmpvars.m
    clear GTC15SpurWMxHL
end
if isDepthSection
    for ii = 1:nFiles
        %Create Transparancy Mask
        WiggleAlpha = sign(RadarDepth{ii}); WiggleAlpha(WiggleAlpha<0) = 0;
        WiggleAlpha = conv2(WiggleAlpha,triang(3),'same')./sum(triang(3));
        WiggleAlpha = WiggleAlpha.*(tukeywin(size(WiggleAlpha,1),.09).*ones(1,size(WiggleAlpha,2)));
        WiggleAlpha(WiggleAlpha<1) = 0;
        % Plot Density, Overlay Peak Amplitudes
        figure();imagesc(Traverse{ii},DepthAxis{ii},1000.*DensityModel{ii});colormap(yet_white);freezeColors;hold on;
        imagesc(Traverse{ii},DepthAxis{ii},sign(RadarDepth{ii}),'AlphaData',WiggleAlpha);colormap([0,0,0]);freezeColors;
        colormap(yet_white);hlay = colorbar; set(hlay,'YDir','reverse','fontsize',14,'fontweight','bold','Ticks',[310,350,400,450,500,550,600,625]);
        set(get(hlay,'ylabel'),'String','Density [kg/m^{3}]', 'rotation', 270,'Units', 'Normalized', 'Position', [4, 0.5, 0])
        title('Density Tomogram')
        xlabel('Distance [km]')
        ylabel('Depth [m]','rotation',270, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0])
        set(gca,'fontsize',14,'fontweight','bold')
        set(gca,'YTick',[0,2.5,5,7.5,10,12.5,15,17.5,20,22.5])
        
        % Plot Density Anamoly, Overlay Peak Amplitudes
        figure();imagesc(Traverse{ii},DepthAxis{ii},1000.*DensityAnomalyModel{ii});colormap(SplitJet);caxis([-15,15]);freezeColors;hold on;
%         plot(Traverse{ii},depth550,'--k','linewidth',3);
        imagesc(Traverse{ii},DepthAxis{ii},sign(RadarDepth{ii}),'AlphaData',WiggleAlpha);colormap([0,0,0]);freezeColors;
        colormap(SplitJet);hlay = colorbar; %set(hlay,'YDir','reverse','fontsize',14,'fontweight','bold');
        set(get(hlay,'ylabel'),'String','Deviation from Mean Density [kg/m^{3}]', 'rotation', 270,'Units', 'Normalized', 'Position', [4, 0.5, 0])
        title('Density Anomaly Tomogram')
        xlabel('Distance [km]')
        ylabel('Depth [m]','rotation',270, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0])
        set(gca,'fontsize',14,'fontweight','bold')
        set(gca,'YTick',[0,2.5,5,7.5,10,12.5,15,17.5,20,22.5])
        
        figure();
        hDA1 = subplot(2,1,1);
        % Plot Mean Average Density Deviation
%         plot(Traverse{ii},1000.*MeanDensityDeviation{ii},'k','linewidth',3)
%         title('Mean Average Deviation - Density [kg/m^{3}]')
        % Plot Surface Density Deviation
        plot(Traverse{ii},1000.*SurfaceDensityDeviation{ii},'k','linewidth',3)
        title('Surface Density - Deviation from Mean [kg/m^{3}]')
        grid on
        set(gca,'fontsize',14,'fontweight','bold')
%         set(hDA1,'units','normalized')
%         hDA1pos = get(hDA1,'Position');
%         set(hDA1,'position',[hDA1pos(1),hDA1pos(2)+.25, hDA1pos(3), hDA1pos(2)-.25])
        
        hDA2 = subplot(2,1,2);
        % Plot Density Anamoly, Overlay Peak Amplitudes
        imagesc(Traverse{ii},DepthAxis{ii},1000.*DensityAnomalyModel{ii});colormap(SplitJet);caxis([-15,15]);freezeColors;hold on;
        imagesc(Traverse{ii},DepthAxis{ii},sign(RadarDepth{ii}),'AlphaData',WiggleAlpha);colormap([0,0,0]);freezeColors;
        colormap(SplitJet);hlay = colorbar; %set(hlay,'YDir','reverse','fontsize',14,'fontweight','bold');
        set(get(hlay,'ylabel'),'String','Deviation from Mean Density [kg/m^{3}]', 'rotation', 270,'Units', 'Normalized', 'Position', [4, 0.5, 0])
        title('Density Anomaly Tomogram')
        xlabel('Distance [km]')
        ylabel('Depth [m]','rotation',270, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0])
        set(gca,'fontsize',14,'fontweight','bold')
        set(gca,'YTick',[0,2.5,5,7.5,10,12.5,15,17.5,20,22.5])
        
        % Plot Depth-Age Tomography
        figure();imagesc(Traverse{ii},DepthAxis{ii},AgeModel{ii});colormap(yet_white);freezeColors;hold on;
%         imagesc(Traverse{ii},DepthAxis{ii},sign(RadarDepth{ii}),'AlphaData',WiggleAlpha);colormap([0,0,0]);freezeColors;
        colormap(yet_white);hlay = colorbar; set(hlay,'YDir','reverse','fontsize',14,'fontweight','bold','Ticks',[Ages{ii}]);
        for kk = 1:size(DepthAge{ii},2)
        plot(Traverse{ii},DepthAge{ii}(:,kk),'k','linewidth',3)
        end
        set(get(hlay,'ylabel'),'String','Age [a]', 'rotation', 270,'Units', 'Normalized', 'Position', [4, 0.5, 0])
        set(gca,'YTick',[0,2.5,5,7.5,10,12.5,15,17.5,20,22.5])
        title('Isochrone Tomogram')
        xlabel('Distance [km]')
        ylabel('Depth [m]','rotation',270, 'Units', 'Normalized', 'Position', [-0.08, 0.5, 0])
        set(gca,'fontsize',14,'fontweight','bold')        
        
    end
    %% Image Depth Section
    for ii = 1:nFiles
        figure();imagesc(Traverse{ii},DepthAxis{ii},AGCgain(RadarDepth{ii},size(RadarDepth{ii},1)./round(3.5),2));
        colormap(Smoke);%hold on;
%         for kk = 1:size(DepthAge{ii},2)
%         plot(Traverse{ii},DepthAge{ii}(:,kk),'k','linewidth',3)
%         end
        title('Core 15 Spur West - Depth Section')
        xlabel('Distance [km]')
        ylabel('Depth [m]','rotation',270, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0])
        set(gca,'fontsize',14,'fontweight','bold')
    end
end
%% Save HVA Output
isSaveHVA = 0;
if isSaveHVA
    % Save Rough Results
    HVAfilename = '6-12-17-Core15-Spur-W-HVA-SurfaceForcing-2016Corrected.mat';
    save(HVAfilename,'AirTo','DirectTo','DirectToVar','deltaT',...
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
    HVAsmoothFilename = '6-12-17-Core15-Spur-W-HVAsmooth-SurfaceForcing-2016Corrected.mat';
    save(HVAsmoothFilename,'dhTWT','dhTWTvar','dhSnowWaterEqv',...
        'dhSnowWaterEqvVar','dhDensity','dhDensityVar','dhDepth',...
        'dhDepthVar','TWT','TWTvar','SnowWaterEqv',...
        'SnowWaterEqvVar','Density','DensityVar','Depth','DepthVar','LayerSnowWaterEqv',...
        'LayerSnowWaterEqvVar','LayerDensity','LayerDensityVar','LayerThickness',...
        'LayerThicknessVar','SurfaceDensity','SurfaceDensityVar',...
        'SurfaceDensification','SurfaceDensificationVar','ForcingDensity',...
        'ForcingDensityVar','ForcingDepth','ForcingDepthVar','AverageAccumulation''-v7.3');
    
    % Save Bootstrap Distributions (Feed into MMxHL Modeling)
    BootstrapFilename = '6-12-17-Core15-Spur-W-Bootstraps-2016Corrected.mat';
    save(BootstrapFilename,'xToRef','xRhoDir','xRhoRef','xDepthDir','xDepth','-v7.3');
end
    % Save Travel Time Picks
    if isWriteTimeHorizons
%         save('6-2-16-Core7-Spur-W-TimeHorizon.mat','DirectFBpick','ReflectionFBpick','-v7.3');
save('6-12-17-Core15-Spur-W1-Surface-8chan-TimeHorizon.mat','DirectFBpick','exportDirectTravelTimes',...
    'ReflectionFBpick','exportReflectionTravelTimes','-v7.3');
    end
    
    % Save Modeled Output
    % Must Include GPS!
    isSaveMxHL = 0;
    if isSaveMxHL
        MxHLFilename = 'GTC15SpurWMxHL.mat';
        GTC15SpurWMxHL = struct('Traverse',{Traverse},'DepthAxis',{DepthAxis},...
            'RadarDepth',{RadarDepth},'AgeModel',{AgeModel},...
            'DensityModel',{DensityModel},'DensityAnomalyModel',{DensityAnomalyModel},...
            'AvgDensityModel',{AvgDensityModel},'MeanDensityDeviation',...
            {MeanDensityDeviation},'SurfaceDensityDeviation',{SurfaceDensityDeviation},...
            'IsoChrones', {Isochrones},'Ages',{Ages},'AverageAccumulation',{AverageAccumulation});
        save(MxHLFilename,'GTC15SpurWMxHL','-v7.3');

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
        axis ij
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
        axis ij
        grid on
        title('Average Density [kg/m^{3}]')
        subplot(3,1,3)
        for dh = 2:nDirectHorizon
            shadedErrorBarT8([],dhDepth{dh,ii},...
                sqrt(dhDepthVar{dh,ii}),1,{'Color',[0,0,0],'linewidth',1.5});
            hold on;
        end
        freezeColors
        axis tight
        axis ij
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
        axis ij
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
        axis ij
        grid on
        title('Average Density [kg/m^{3}]')
        subplot(3,1,3)
        for rh = 1:nReflectionHorizon
            shadedErrorBarT8([],Depth{rh,ii},...
                sqrt(DepthVar{rh,ii}),1,{'Color',[0,0,0],'linewidth',1.5});
            hold on;
        end
        freezeColors
        axis tight
        axis ij
        grid on
        title('Depth [m]')
        set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,...
            'FontWeight','Bold', 'LineWidth', 1);
    end
end
%% Joint Figure for Snow Surface Data
if isSWEDISH
    % plot Direct Wave Data
    for ii = 1:nFiles
        distance = Traverse{ii};
        figure();
        subplot(3,1,1)
%         for dh = 2:nDirectHorizon
%             shadedErrorBarT8(distance,dhSnowWaterEqv{dh,ii},...
%                 sqrt(dhSnowWaterEqvVar{dh,ii}),1,{'Color',[0,0,0],'linewidth',1.5});
%             hold on;
%         end
        for rh = 1:nReflectionHorizon
            shadedErrorBarT8(distance,AverageAccumulation{ii},...
                sqrt(SnowWaterEqvVar{rh,ii}),1,{'Color',[0,0,0],'linewidth',1.5});            
%             shadedErrorBarT8(distance,SnowWaterEqv{rh,ii},...
%                 sqrt(SnowWaterEqvVar{rh,ii}),1,{'Color',[0,0,0],'linewidth',1.5});
            hold on;
        end
        freezeColors
        axis tight
        axis ij
        grid on
        title('Accumulation [m w.e.]')
        subplot(3,1,2)
        for dh = 2:nDirectHorizon
            shadedErrorBarT8(distance,dhDensity{dh,ii}.*1000,...
                sqrt(dhDensityVar{dh,ii}).*1000,1,{'Color',[0.5,0,0],'linewidth',1.5});
            hold on;
        end
        for rh = 1:nReflectionHorizon
            shadedErrorBarT8(distance,Density{rh,ii}.*1000,...
                sqrt(DensityVar{rh,ii}).*1000,1,{'Color',[0.5,0,0],'linewidth',1.5});
            hold on;
        end
        freezeColors
        axis ij
        axis tight
        grid on
        title('Average Density [kg/m^{3}]')
        subplot(3,1,3)
%         for dh = 2:nDirectHorizon
%             shadedErrorBarT8(distance,dhDepth{dh,ii},...
%                 sqrt(dhDepthVar{dh,ii}),1,{'Color',[1,0.81,0],'linewidth',1.5});
%             hold on;
%         end
        for rh = 1:nReflectionHorizon
            shadedErrorBarT8(distance,Depth{rh,ii},...
                sqrt(DepthVar{rh,ii}),1,{'Color',[1,0.81,0],'linewidth',1.5});
            hold on;
        end
        freezeColors
        axis ij
        axis tight
        grid on
        title('Depth [m]')
        xlabel('Distance [km]')
        set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,...
            'FontWeight','Bold', 'LineWidth', 1);
    end
end