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

% GreenTrACS 2017
% dataDir = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/rawNC/6-18-17-SummitRouteI';
% dataDir = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/rawNC/6-17-17-Core16Spiral';
% dataDir = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/rawNC/6-16-17-Core16SpurW';
% dataDir = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/rawNC/6-15-17-Core15Core16Traverse';
% dataDir = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/rawNC/6-14-17-Core15Spiral';
dataDir = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/rawNC/6-13-17-Core15SpurW';
% dataDir = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/rawNC/6-12-17-Core14Core15Traverse';
% dataDir = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/rawNC/6-11-17-Core14Spiral';
% dataDir = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/rawNC/6-10-17-Core14SpurW';
% dataDir = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/rawNC/6-9-17-Core13Core14Traverse';
% dataDir = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/rawNC/6-8-17-Core13SpurE';
% dataDir = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/rawNC/6-6-17-Core12Core13Traverse';
% dataDir = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/rawNC/6-3-17-Core12SpurE';
% dataDir = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/rawNC/6-2-17-Core12SpurW';
% dataDir = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/rawNC/6-1-17-Core11Core12Traverse';
% dataDir = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/rawNC/5-31-17-Core11SpurW';
% dataDir = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/ArcticDataCenter/PulseEKKO/500MHz/rawNC/5-27-17-Core11SpurE';

% Add additional useful pathways
addpath '/home/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/Save'
addpath '/home/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/Save/matData'
addpath '/home/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/TM'
addpath '/home/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/time2depth'
addpath './functions';
addpath './supplementalData';
workDir = pwd;

%% Control Floe
% Parallel Computing Enabled
isParallel = 1;

% Read Data
isReadNC = 1;                  % Read Multiplexed Data
isLoadTimeHorizons = 1;        % Load Previously Picked Time Horizons
isPickTravelTimeHorizons = 0;  % Pick Travel-Time Horizons
isLoadIRH = 1;                 % Load Previously Picked IRHs
isPickAgeHorizons = 0;         % Pick Age Horizons
isLoadDepthHorizons = 1;       % Load Previously Picked Depth Horizons
isPickDepthHorizons = 0;       % Pick Isochronous Depth Horizons
isLoadHVA = 1;                 % Load Previous Horizon Velocity Analysis
isLoadMxHL = 0;                % Load Previous MxHL Model Results
isLoadGPS = 1;                 % Load GPS for MxRadar
isGreenTracsFirnCore = 1;      % Load GreenTracs Firn Core Data
isMEaSUREs = 0;                % Load NASA MEaSUREs Surface Velocity

% Export Data
isWriteTimeHorizons = 0;% Save Travel-Time Picks
isSaveHVA = 0;          % Save Horizon Velocity Analysis
isSaveMxHL = 1;         % Save Modeled Output

% Process Data
isReduceData = 1;       % Remove Every nth Trace from Data Gather
isTrimTWT = 0;          % Truncate Recorded Data for Near-Surface Analysis
isKill = 1;             % Kill Unanted Channels
isMedianSubtraction = 1;% Background Median Subtraction Filter
isSWEDISH = 1;          % Perform Surface Velocity Analysis
isDepthSection = 1;     % NMO Correction, Stacking, Depth Conversion, Image
isFXdecon = 1;          % Fx-Predicitive Deconvolution (Depth Domain)

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
yet_black = [[1,1,1];[.9463,.9463,1];yet_white(2,:);[(yet_white(2,1:2)-yet_white(3,1:2))./2 + yet_white(3,1:2),...
    1];yet_white(3:end,:);[linspace(yet_white(end,1)-.0625,0,3)',zeros(3,1),zeros(3,1)]];

    % Load ColorBrewer 
    set(0,'DefaultAxesFontName','FreeSerif')
    set(0,'DefaultTextFontName','FreeSerif')
    [colorbrew , ~, ~] = brewermap(256,'RdYlBu');
    colorbrew=flipud(colorbrew);
    colorQ = colorbrew(round(quantile(1:256,[0.85,0.5,0.15])),:);
    c1 = colorQ(1,:);
    c2 = colorQ(2,:);
    c3 = colorQ(3,:);
    

TraverseDistance = [15,15,15];  % Approx. Distance of Radar Files [km]
fileNames = dir([dataDir,'/','*.nc']);
lineNo = [0,1,2,4,7];%[3,4,5];               % Array of data "LINE" numbers
nFiles = length(lineNo);        % Number of Files
nChan = 9;                      % Number of Recorded Channels
chan =  1:nChan;                % Linear Array of Record Channels
liveChan = chan;

% Establish Tx Rx Geometry for CMP Gathering
nTx = 3;               % Number of Transmitters in Sequence
nRx = 3;               % Number of Receivers in Sequence

% Allocate Memory
time = cell(1,nFiles);

%% Load HVA Results
if isLoadHVA
    cd '/home/tatemeehan/GreenTracs2017/MXHL';
    load('Core15SpurW_HVAsurfaceForcing_061919.mat');
    load('Core15SpurWHVA_061919.mat');
    load('Core15SpurW_Bootstraps_061919.mat');
    cd(workDir)
end

%% Import GPS Information
addpath '/home/tatemeehan/GreenTracs2017/GPS/Core15SpurW061317';
if isLoadGPS
    load('MxRadarGPSCore15SpurW061317.mat')
end
%% Read GPR Data
if isReadNC

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

end

%% Branch Processing for Direct and Reflected Arrivals
directRadar = Radar;
%% Snow Water Equivalent and Density for Ice Sheet Height

    if isLoadTimeHorizons

        TimeHorizonFilename = 'Core15SpurWCompleteTimeHorizon.mat';

        load(TimeHorizonFilename);
        
        % Determine Number of Direct Wave Horizons
        nDirectHorizon = size(DirectFBpick,2);
        % Determine Number of Reflection Horizons
        nReflectionHorizon = size(ReflectionFBpick,2);
   
        % Plot Common-offset gathers
        isPlotOffsetGathers = 0;
        isTopMute = 0;
        
        VisualizeDataGathers
        
    end
    
%% Semi-Automatic Radar Wave Picking
if isSWEDISH && ~isLoadHVA
    if isPickTravelTimeHorizons
    fprintf('Begin PolarPicker \n')
    display(' ')
    tic

    PickTravelTimeHorizons    

    fprintf('PolarPicker Done \n')
    toc
    display(' ')
    end
%% Perform Horizon Velocity Analysis for Estimates of Density, Depth, & SWE
    fprintf('Begin Horizon Velocity Analysis \n')
    tic

    % Bias Calibration at Core 15W Pit (Observed - Estimated);
    % This is a 1.72% Adjustment
    % Correction has been removed
    velocityBias = 0.004;      % [m/ns]
    densityBias = -.0258902;   % [g/cm3]
    
    if isGreenTracsFirnCore
        isCoreDepthAge = 1; % Use Age Depth Profile from Local Core Site
        coreNo = [15]; % Array of Firn Cores 1-16 to include in analysis
        depthAgeFilename = 'Core15_age_scale.txt';
        depthDensityFilename = 'Core15_depth_density.txt';
        
        % Annual Accumulation Correction from Firn Core Chemistry [mwe]
        GreenTracsFirnCore
    else
        % Assume average is 0.3 [mwe]
        coreAccumulation = 0.3;
        % Average over 2.5 years Winter 2015 - Summer 2017
        ageInterval = 2.5;
        % If Radar is not Paired with GPS or Ice Core data is not used 
        % assume Ice Core is located adjacent to the first k traces of the 
        % first radar file
        k = 100; % Number of Radar Estimates to Estimate Core Site Average
        iceCoreFileIx = 1;
        iceCoreIx = 1:k;
        warning('Undefined Firn Core!.. Using default accumulation 0.3 mwe and age 2.5 yrs.')
    end
    
    % Toggle Inversion Scheme
    % Air-Coupled Wave Solution
    isL1Air = 0;
    isL2Air = 1;
    % Surface-Coupled Wave
    isL1LMO = 0;
    isL2LMO = 1;
    % Reflected Wave
    isL1NMO = 0;
    isL2NMO = 1;
    
    % Memory Allocation
    HVAmemoryAllocation
    
    % Run Horizon Velocity Analysis
    HorizonVelocityAnalysis
elseif isGreenTracsFirnCore
    isCoreDepthAge = 1; % Use Age Depth Profile from Local Core Site
    coreNo = [15]; % Array of Firn Cores 1-16 to include in analysis
    depthAgeFilename = 'Core15_age_scale.txt';
    depthDensityFilename = 'Core15_depth_density.txt';
    
    % Annual Accumulation Correction from Firn Core Chemistry [mwe]
    GreenTracsFirnCore
end   
%% Apply Residual Trace Shifts to Data Gathers
if isSWEDISH || isLoadHVA && ~isLoadMxHL

    DataStaticShift
    
    fprintf('Horizon Velocity Analysis Done \n')
    toc
    display(' ')

end
%% Extrapolate Surface Density from Radar Forcing

    MxHLforcingParameter

%% Create Depth Section and Density Model
if isDepthSection && ~isLoadMxHL
    %% Stacking Velocity Model
    fprintf('Begin Stacking Velocity Extrapolation \n')
    tic
    
    VelocityModel
    
    fprintf('Velocity & Density Extrapolation Done \n')
    toc
    display(' ')
    

    %% Normal Moveout Correction
    fprintf('Begin Normal Moveout Correction \n')
    tic
    
    NormalMoveoutCorrection    
    
    fprintf('Normal Moveout Correction Done \n')
    toc
    display(' ')
    %% Stack Gathers
    fprintf('Begin Offset-Gather Stacking \n')
    tic

    StackGathers
    
    fprintf('Stacking Done \n')
    toc
    display(' ')
    %% Post-Stack FX-Deconvolution
    isPostStackFXdecon = 0;
    if isPostStackFXdecon
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
    %% Time to Depth Conversion
    isTime2Depth = 0;
    if isTime2Depth
    fprintf('Begin Time-Depth Conversion \n')
    tic
    
    TimeToDepthConversion
    
    fprintf('Time-Depth Conversion Done \n')
    toc
    display(' ')
    end
    %% Resize Density Model
    fprintf('Begin Regrid of MxHL Model \n')
    tic

    MxHLmodel

    fprintf('MxHL Model Regrid Done \n')
    toc
    display(' ')
    %% Radar Time to Stratigraphic Age Image
    fprintf('Begin Time-Age Conversion \n')
    tic
    
    TimeToAgeConversion

    fprintf('Time-Age Conversion Done \n')
    toc
    display(' ')
    
    %% PickAge Horizons for Residual Update
    if isPickAgeHorizons
        fprintf('Begin Picking Age Image \n')
        tic
        
        PickAgeHorizons
        
        % Write Isochrone Picks to .mat
        isWriteIRH =0;
        if isWriteIRH
            cd '/home/tatemeehan/GreenTracs2017/MXHL';
            save('isochronePicksCore15SpurW110319.mat','isochronePick','-v7.3')
            cd(workDir)
        end
        
        fprintf('Depth-Deposition Picking Done \n')
        toc
        display(' ')
        
    elseif isLoadIRH
        cd '/home/tatemeehan/GreenTracs2017/MXHL';
        IRH = load('isochronePicksCore15SpurW070719.mat');
%         IRH = load('isochronePicksCore15SpurW110319.mat');
        isochronePick = IRH.isochronePick;
        clear IRH
        cd(workDir)
        fprintf('Isochrone Reflection Horizon Picks Loaded \n')
        disp(' ')
    end
    
    %% Isochrone Model Update and Trace Flattening
    if isPickAgeHorizons || isLoadIRH
        fprintf('Begin Isochrone Model Update & Trace Flattening \n')
        tic
        % Calculate Stratigraphic Age Residual
        CalculateAgeResidual
        % Update Model with Perturbations
        UpdateAgeHorizons
        
        fprintf('Isochrone Model Update & Trace Flattening Done \n')
        toc
        display(' ')
    end
    %% Deposition Image FX-Deconvolution
    if isFXdecon
        fprintf('Begin Deposition Image FX-Predictive Deconvolution \n')
        tic
        
        % Select Depth or Spatial Deconvolution Gates: 1 on; 0 off;
        isDepthFXdecon = 1;
        isSpatialFXdecon = 0;
        
        DepositionFXdeconTX
        
        fprintf('Deposition Image FX-Deconvolution Done \n')
        toc
        display(' ')
    end
    
    %% Isochrone Model Update and Trace Flattening
%     if isPickAgeHorizons || isLoadIRH
%         fprintf('Begin Trace Unflattening \n')
%         tic
%         
%         DepositionUnflattening
%         
%         fprintf('Trace Unflattening Done \n')
%         toc
%         display(' ')
%     end
    
    %% Radar Stratigraphic Age to Depth Image
    fprintf('Begin Age-Depth Conversion \n')
    tic
      AgeToDepthConversion

    fprintf('Age-Depth Image Done \n')
    toc
    display(' ')
end
%% Pick Depth Horizons for Residual Update
if isPickDepthHorizons
    fprintf('Begin Picking Depth Image \n')
    tic
    
    PickDepthHorizons
    
    % Write Isochrone Picks to .mat
    isWriteDepthHorizons = 0;
    if isWriteDepthHorizons
        cd '/home/tatemeehan/GreenTracs2017/MXHL';
        save('depthPicksCore15SpurMaster072919.mat','depthPick','-v7.3')
        cd(workDir)
    end
    
    fprintf('Depth Domain Isochrone Picking Done \n')
    toc
    display(' ')
    
elseif isLoadDepthHorizons
    cd '/home/tatemeehan/GreenTracs2017/MXHL';
    IRH = load('depthPicksCore15SpurMaster072919.mat');
    depthPick = [IRH.depthPick];
%     clear IRH
    cd(workDir)
    fprintf('Isochrone Reflection Horizon Picks Loaded \n')
    disp(' ')
end
    
    %% Isochrone Model Update and Trace Flattening
%     if isPickDepthHorizons || isLoadDepthHorizons
%         fprintf('Begin Isochrone Model Update & Trace Flattening \n')
%         tic
%         % Calculate Stratigraphic Age Residual
        ReCalculateAgeResidual        
%         % Update Model with Perturbations
%         UpdateAgeHorizons
%         
%         fprintf('Isochrone Model Update & Trace Flattening Done \n')
%         toc
%         display(' ')
%     end

%% Save HVA Output
if isSaveHVA
    % Save Rough Results
    HVAfilename = 'Core15SpurWHVA_070819.mat';
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
    HVAsmoothFilename = 'Core15SpurW_HVAsurfaceForcing_061919.mat';
    save(HVAsmoothFilename,'dhTWT','dhTWTvar','dhSnowWaterEqv',...
        'dhSnowWaterEqvVar','dhDensity','dhDensityVar','dhDepth',...
        'dhDepthVar','TWT','TWTvar','SnowWaterEqv',...
        'SnowWaterEqvVar','Density','DensityVar','Depth','DepthVar','LayerSnowWaterEqv',...
        'LayerSnowWaterEqvVar','LayerDensity','LayerDensityVar','LayerThickness',...
        'LayerThicknessVar','SurfaceDensity','SurfaceDensityVar',...
        'SurfaceDensification','SurfaceDensificationVar','ForcingDensity',...
        'ForcingDensityVar','ForcingDepth','ForcingDepthVar','AverageAccumulation','-v7.3');
    
    % Save Bootstrap Distributions (Feed into MMxHL Modeling)
    BootstrapFilename = 'Core15SpurW_Bootstraps_061919.mat';
    save(BootstrapFilename,'xToRef','xRhoDir','xRhoRef','xDepthDir','xDepth','-v7.3');
end
    % Save Travel Time Picks
    if isWriteTimeHorizons
%         save('6-2-16-Core7-Spur-W-TimeHorizon.mat','DirectFBpick','ReflectionFBpick','-v7.3');
%'6-12-17-Core15-Spur-W1-Surface-8chan-TimeHorizon.mat'
cd '/home/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/Save/matData'
save('Core15SpurW080219CompleteTimeHorizon.mat','DirectFBpick','exportDirectTravelTimes',...
    'ReflectionFBpick','exportReflectionTravelTimes','-v7.3');
cd(workDir)

    end
    
    % Save Modeled Output
    if isSaveMxHL
        MxHLFilename = 'GTC15SpurWMxHL_042820.mat';
        GTC15SpurWMxHL = struct('DistanceAxis',{Traverse},'DepthAxis',{DepthAxis},...
            'RadarDepth',{RadarDepth},'TimeAxis',{TimeAxis},'RadarStack',{RadarStack},'DepositionAxis',{DepositionAxis},'RadarDeposition',{RadarDeposition},'AgeModel',{bestAgeModel},...
            'DensityModel',{DensityModel},'DensityAnomalyModel',{DensityAnomalyModel},...
            'AvgDensityModel',{AvgDensityModel},'MeanDensityDeviation',...
            {MeanDensityDeviation},'SurfaceDensityDeviation',{SurfaceDensityDeviation},...
            'IsochroneIx',{IRH.depthPick},'DepthIsoChrones', {depthPick},'AgeIsochrones',{agePick},'SMBmodel',{instantSMB},'AverageAccumulation',{AverageAccumulation3},...
            'VarAccumulation',{varAccumulation3},'iceCoreIx',{iceCoreIx});
        cd '/home/tatemeehan/GreenTracs2017/MXHL/'
        save(MxHLFilename,'-struct','GTC15SpurWMxHL','-v7.3');
        cd(workDir)
    end
    
%% Radar Imagery and Snow Data Visualization
MakeFigures