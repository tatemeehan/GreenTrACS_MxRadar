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
addpath './functions';
addpath './supplementalData';
%% Control Floe
% Parallel Computing Enabled
isParallel = 1;

% Read Data
isReadSensorsSoftware = 1;     % Read Multiplexed Data
isLoadTimeHorizons = 1;        % Load Previously Picked Time Horizons
isPolarPicker = 0;             % Pick Travel-Time Horizons
isLoadHVA = 1;                 % Load Previous Horizon Velocity Analysis
isLoadMxHL = 0;                % Load Previous MxHL Model Results
isLoadGPS = 1;                 % Load GPS for MxRadar

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
fileNames = dir([dataDir,'/','*.nc']);
lineNo = [0,1,2];%,4,7];               % Array of data "LINE" numbers
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
    load('6-12-17-Core15-Spur-W-HVA-SurfaceForcing-Corrected.mat');
    load('6-12-17-Core15-Spur-W-HVAsmooth-SurfaceForcing-Corrected.mat');
    load('6-12-17-Core15-Spur-W-Bootstraps-Corrected.mat');
end

%% Import GPS Information
addpath '/home/tatemeehan/GreenTracs2017/GPS/Core15SpurW061317';
if isLoadGPS
    load('MxRadarGPSCore15SpurW061317.mat')
end
%% Read GPR Data
if isReadSensorsSoftware

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

        TimeHorizonFilename = '6-12-17-Core15-Spur-W1-Surface-8chan-TimeHorizon.mat';
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
if isSWEDISH
    if isPolarPicker
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
    
    HorizonVelocityAnalysis

end   
%% Apply Residual Trace Shifts to Data Gathers
if isSWEDISH || isLoadHVA && ~isLoadMxHL

    DataStaticShift
    
    fprintf('Horizon Velocity Analysis Done \n')
    toc
    display(' ')

end
%% Extrapolate Surface Density from Radar Forcing
% Bias Calibration at Core 15W Pit (Observed - Estimated);
% Applied Within HVA

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
    %% Time to Depth Conversion
    fprintf('Begin Time-Depth Conversion \n')
    tic
    
    TimeToDepthConversion
    
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

    MxHLmodel

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
%         figure();imagesc(Traverse{ii},DepthAxis{ii},RadarDepth{ii});
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