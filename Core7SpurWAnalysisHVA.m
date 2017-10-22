%% Core 7 Analysis
% Tate Meehan GreenTrACS2016
% clear; close all; clc;
addpath '/home/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/TM'
load('6-2-16-Core7-Spur-W-HVA.mat');%load('Core7FieldRecordGuideR0.mat');

%% Core 7 SnoPit Data
    tare = 775/1000;                         % Tare Weight of Kelly Cutter
    % Measurement Intervals
    depthCut = [.10,.20,.30,.40,.50,.60,.71] - 0.05;
    % Snow Weight Measurements
    cut1 = [1070,1165,1154,1145,1180,1230,1221]/1000;
    cut2 = [1068,1140,1131,1160,1150,1180,1211]/1000;
    cut = vertcat(cut1,cut2);
    weights = (cut-tare);          % Tare Subrtraction - Yeilds Density
    Kelly = weights(1,:)-weights(2,:);  % Diff Measurements for Error Est.
    fightingWeight = mean(weights);     % Average Density at Depth
    KellyWeight = mean(fightingWeight); % Bulk Density

%% Input Core 7 Density Log
% Format : Profile = [ Depth, Thickness, Density, bulkDesnity, cumSWE, SWE]
Profile = zeros(30,6);
    % Input Depth [m] and Density [%]
    Profile(1,1) = depthCut(1); Profile(1,3) = fightingWeight(1);
    Profile(2,1) = depthCut(2); Profile(2,3) = fightingWeight(2);
    Profile(3,1) = depthCut(3); Profile(3,3) = fightingWeight(3);
    Profile(4,1) = depthCut(4); Profile(4,3) = fightingWeight(4);
    Profile(5,1) = depthCut(5); Profile(5,3) = fightingWeight(5);
    Profile(6,1) = depthCut(6); Profile(6,3) = fightingWeight(6);
    Profile(7,1) = depthCut(7); Profile(7,3) = fightingWeight(7);
%     Profile(1,1) = max(depthCut); Profile(1,3) = KellyWeight;
    Profile(8,1) = 1.950;    Profile(8,3) = 0.403;
    Profile(9,1) = 2.465;    Profile(9,3) = 0.401;
    Profile(10,1) = 3.600;    Profile(10,3) = 0.438;
    Profile(11,1) = 4.510;    Profile(11,3) = 0.488;
    Profile(12,1) = 5.395;    Profile(12,3) = 0.465;
    Profile(13,1) = 6.425;    Profile(13,3) = 0.456;
    Profile(14,1) = 7.515;    Profile(14,3) = 0.492;
    Profile(15,1) = 8.440;    Profile(15,3) = 0.516;
    Profile(16,1) = 9.445;   Profile(16,3) = 0.527;
    Profile(17,1) = 10.400;  Profile(17,3) = 0.538;
    Profile(18,1) = 11.445;  Profile(18,3) = 0.545;
    Profile(19,1) = 12.195;  Profile(19,3) = 0.568;
    Profile(20,1) = 13.240;  Profile(20,3) = 0.569;
    Profile(21,1) = 14.055;  Profile(21,3) = 0.565;
    Profile(22,1) = 14.990;  Profile(22,3) = 0.587;
    Profile(23,1) = 15.875;  Profile(23,3) = 0.577;
    Profile(24,1) = 16.580;  Profile(24,3) = 0.582;
    Profile(25,1) = 17.425;  Profile(25,3) = 0.583;
    Profile(26,1) = 18.395;  Profile(26,3) = 0.590;
    Profile(27,1) = 19.295;  Profile(27,3) = 0.610;
    Profile(28,1) = 20.170;  Profile(28,3) = 0.615;
    Profile(29,1) = 21.060;  Profile(29,3) = 0.604;
    Profile(30,1) = 21.875;  Profile(30,3) = 0.615;
    
    
    % Compute Bulk Density
    cumDensity = cumsum(Profile(:,3));
    for ii = 1:length(Profile)
    Profile(ii,4) = cumDensity(ii)/ii;  % Bulk Density
    end
    
    % Input Layer Thickness [m]
    Profile(1,2) = Profile(1,1);
    Profile(2:end,2) = diff(Profile(:,1));
    
    % Input SWE Estimate [m]
    for ii = 1:length(Profile)
        Profile(ii,5) = Profile(ii,1).*Profile(ii,4); % Cumulative SWE
        Profile(ii,6) = Profile(ii,2).*Profile(ii,3); % Layer SWE
    end
    
    % Estimate Velocity From Core Log
    for jj = 1:length(Profile)
       % Interval Velocity
    Vint(jj,1) =  DryCrimVRMS( Profile(jj,3) );
    VRMS(jj,1) = DryCrimVRMS( Profile(jj,4) );
    end
    Vint2 = Vint.^2;
    % Estimate Travel Time Through Layer
    for jj = 1:length(Profile)
        ti(jj,1) = (Vint(jj,1)./Profile(jj,2)).^(-1);
    end
    
    % Estimate Average Velocity and RMS Velocity
    for jj = 1:length(Profile)
        Vavg(jj,1) = sum(Profile(1:jj,2))./(sum(ti(1:jj,1)));
%         Vrms(jj,1) = sqrt(sum(Vint(1:jj,1).*(Profile(1:jj,2)))./(sum(ti(1:jj,1))));
        Vrms(jj,1) = sqrt(sum(Vint2(1:jj,1).*ti(1:jj,1))./(sum(ti(1:jj,1))));
    end
    
    % Estimate Interval Velocity
    toAVG = 2.*(Profile(:,1)./Vavg);
    toRMS = 2.*(Profile(:,1)./Vrms);

    [ VintAvg, ~, ZintAvg ] = Dix( Vavg, toAVG );
    [ VintRMS, ~, ZintRMS ] = Dix( Vavg, toRMS );
    % Error Analysis
    ErrVintAvg = sqrt(mean((Vint - VintAvg').^2));
    ErrVintRMS = sqrt(mean((Vint - VintRMS').^2));
    ErrVavgVrms = sqrt(mean((Vavg - Vrms).^2));
    ErrVRMS = sqrt(mean((VRMS - Vrms).^2));
    ErrVavg = sqrt(mean((VRMS - Vavg).^2));
    
    %% Core Spur W HVA Data at Core 7
    % Surface Waves
    LMOto = cell2mat(DirectTo);
    LMOtoVar = cell2mat(DirectToVar);    
    VlmoHVA = cell2mat(DirectVelocity);
    VlmoHVAvar = cell2mat(DirectVelocityVar);    
    ZlmoHVA = cell2mat(DirectDepth);
    ZlmoHVAvar = cell2mat(DirectDepthVar);    
    RHOlmoHVA = cell2mat(DirectDensity);
    RHOlmoHVAvar = cell2mat(DirectDensityVar);    
    SWElmo = cell2mat(DirectSWE);
    SWElmoVar = cell2mat(DirectSWEvar);
    % Reflected Waves
    NMOto = cell2mat(ReflectionTo);
    NMOtoVar = cell2mat(ReflectionToVar);    
    VrmsHVA = cell2mat(ReflectionVelocity);
    VrmsHVAvar = cell2mat(ReflectionVelocityVar);    
    VintHVA = cell2mat(IntervalVelocity);    
    VintHVAvar = cell2mat(IntervalVelocityVar);
    ZnmoHVA = cell2mat(ReflectionDepth);
    ZnmoHVAvar = cell2mat(ReflectionDepthVar);
    HintHVA = cell2mat(IntervalThickness);
    HintHVAvar = cell2mat(IntervalThicknessVar);
    RHOrms = cell2mat(ReflectionDensity);
    RHOrmsVar = cell2mat(ReflectionDensityVar);
    RHOint = cell2mat(IntervalDensity);
    RHOintVar = cell2mat(IntervalDensityVar);
    SWErms = cell2mat(SWE);
    SWErmsVar = cell2mat(SWEvar);
    SWEh = cell2mat(SWEint);
    SWEhVar = cell2mat(SWEintVar);
    
    % Grab Near Core 7 Estimates
    nearIx = 1:length(NMOto);  % Shot Gather Ix for Near Core Site Estimates
    % Depth
    RadarDepth = mean([abs(ZlmoHVA(2,nearIx));ZnmoHVA(:,nearIx)],2);
    RadarDepthVar = mean([abs(ZlmoHVAvar(2,nearIx));ZnmoHVAvar(:,nearIx)],2);
    
    for jj = 1:length(RadarDepth)
        [~,CoreValidIx(jj)] = min(abs(Profile(:,1)-RadarDepth(jj)));
    end
    
    CoreDepth = Profile(CoreValidIx,1);
    
    % Cumulative Average Density
    RadarRMSdensity = mean([abs(RHOlmoHVA(2,nearIx));RHOrms(:,nearIx)],2);
    RadarRMSdensityVar = mean([abs(RHOlmoHVAvar(2,nearIx));RHOrmsVar(:,nearIx)],2);
    CoreRMSdensity = Profile(CoreValidIx,4);
    
    % Interval Density
    RadarIntervalDensity = nanmean([abs(RHOlmoHVA(2,nearIx));RHOint(:,nearIx)],2);
    RadarIntervalDensityVar = nanmean([abs(RHOlmoHVAvar(2,nearIx));RHOintVar(:,nearIx)],2);
%     CoreIntervalDensity = Profile(CoreValidIx,3);
    
    % Average Core Intervals to match resolution of radar
    CoreIntervalIx = find(diff(CoreValidIx)>1)+1;
    CoreIntervalTop = CoreValidIx([CoreIntervalIx-1])+1;
    CoreIntervalBot = CoreValidIx(CoreIntervalIx);
    intIx = 0;
    for jj = 1:length(CoreValidIx)
        if ~ismember(jj,CoreIntervalIx)
            CoreIntervalDensity(jj) = Profile(CoreValidIx(jj),3);
        else
            intIx = intIx+1;
            CoreIntervalDensity(jj) = mean(Profile(CoreIntervalTop(intIx):CoreIntervalBot(intIx),3));
        end
    end
    
    % Interval Thickness
    RadarIntervalThickness = nanmean([abs(ZlmoHVA(2,nearIx));HintHVA(:,nearIx)],2);
    RadarIntervalThicknessVar = nanmean([abs(ZlmoHVAvar(2,nearIx));HintHVAvar(:,nearIx)],2);
    % Average Core Intervals to match resolution of radar
    intIx = 0;
    for jj = 1:length(CoreValidIx)
        if ~ismember(jj,CoreIntervalIx)
            CoreIntervalThickness(jj) = Profile(CoreValidIx(jj),2);
        else
            intIx = intIx+1;
            CoreIntervalThickness(jj) = sum(Profile(CoreIntervalTop(intIx):CoreIntervalBot(intIx),2));
        end
    end
%     CoreIntervalThickness = [Profile(CoreValidIx(1),2);diff(Profile(CoreValidIx,1))];
    
    % Average Accumulation
    RadarRMSaccumulation = mean([abs(SWElmo(2,nearIx));SWErms(:,nearIx)],2);
    RadarRMSaccumulationVar = mean([abs(SWElmoVar(2,nearIx));SWErmsVar(:,nearIx)],2);
    CoreRMSaccumulation = Profile(CoreValidIx,5);
    % Interval Accumulation
    RadarIntervalAccumulation = nanmean([abs(SWElmo(2,nearIx));SWEh(:,nearIx)],2);
    RadarIntervalAccumulationVar = nanmean([abs(SWElmoVar(2,nearIx));SWEhVar(:,nearIx)],2);
    RadarDiffAccumulation = [RadarRMSaccumulation(1);diff(RadarRMSaccumulation)];
        % Average Core Intervals to match resolution of radar
    intIx = 0;
    for jj = 1:length(CoreValidIx)
        if ~ismember(jj,CoreIntervalIx)
            CoreIntervalAccumulation(jj) = Profile(CoreValidIx(jj),6);
        else
            intIx = intIx+1;
            CoreIntervalAccumulation(jj) = sum(Profile(CoreIntervalTop(intIx):CoreIntervalBot(intIx),6));
        end
    end
    
    %% Create Scatter Plots
% red1 = [0.85,0.33,0.1];
% orange1 = [1, 0.6,0];
% yellow1 = [1,0.84,0];
% green2 = [0,0.8,0];
% green1 = [0,.7,.4];%[0.23,0.44,0.34];
% blue2 = [0,0.8,0.8];
% blue1 = [0.31,0.4,0.58];
% violet2 = [0.6,0.6,1];
% violet1 = [0.49,0.18,0.56];
% maj3 = [0.8,0.6,1];
% maj2 = [1,.6,.8];
% maj1 = [1,.2,1];
% colors = [red1;orange1;yellow1;green2;green1;blue2;blue1;violet2;violet1;maj3;maj2;maj1];
cmap = colormap('jet');
% modcolors = mod(length(cmap),length(CoreValidIx));
% ncolors = (length(cmap)-modcolors)./length(CoreValidIx);
ncolors = length(CoreValidIx);
colorIx = round(linspace(1,length(cmap),ncolors));
colors = cmap(colorIx,:);
% if length(cmap) > 64
%     colors = cmap(1:ncolors:(end-modcolors),:);
% else
%     colors = cmap(1:ncolors:(end-modcolors),:);
% end
colors(8,1) = colors(8,1).*.95;
figure();
% Depth
subplot(3,2,1);
for jj = 1:length(RadarDepth)
    hold on;
    plot(CoreDepth(jj),RadarDepth(jj),'o','markerfacecolor',colors(jj,:),'markeredgecolor',[0.75,0.75,0.75],'linewidth',.5,'markersize',7)
    % scatter(CoreDepth,RadarDepth,50,colors,'filled','markeredgecolor',[0.75,0.75,0.75],'linewidth',.5);
end
legend('Horizon 1','Horizon 2','Horizon 3','Horizon 4','Horizon 5','Horizon 6','Horizon 7','Horizon 8','Horizon 9','Horizon 10','Horizon 11','Horizon 12','location','westoutside')

plot([0,20],[0,20],'k','linewidth',2); hold on;errorbar(CoreDepth,RadarDepth,sqrt(RadarDepthVar),'ok','markersize',0.5,'linewidth',2);
for jj = 1:length(RadarDepth)
    hold on;
    plot(CoreDepth(jj),RadarDepth(jj),'o','markerfacecolor',colors(jj,:),'markeredgecolor',[0.75,0.75,0.75],'linewidth',.5,'markersize',7)
    % scatter(CoreDepth,RadarDepth,50,colors,'filled','markeredgecolor',[0.75,0.75,0.75],'linewidth',.5);
end
title('Depth [m]')
set(gca,'fontsize',12,'fontweight','bold')
axis tight
grid minor
% RMS Average Density
subplot(3,2,3);
plot([250,550],[250,550],'k','linewidth',2); hold on;errorbar(1000.*CoreRMSdensity,1000.*RadarRMSdensity,1000.*sqrt(RadarRMSdensityVar),'ok','markersize',0.5,'linewidth',2);...
scatter(1000.*CoreRMSdensity,1000.*RadarRMSdensity,50,colors,'filled','markeredgecolor',[0.75,0.75,0.75],'linewidth',.5);
title('Average Density [kgm^{-3}]')
set(gca,'fontsize',12,'fontweight','bold')
axis tight
grid minor
% Interval Thickness
subplot(3,2,2);
plot([0,5],[0,5],'k','linewidth',2); hold on;errorbar(CoreIntervalThickness,RadarIntervalThickness,sqrt(RadarIntervalThicknessVar),'ok','markersize',0.5,'linewidth',2);...
scatter(CoreIntervalThickness,RadarIntervalThickness,50,colors,'filled','markeredgecolor',[0.75,0.75,0.75],'linewidth',.5);
title('Interval Thickness [m]')
set(gca,'fontsize',12,'fontweight','bold')
axis tight
grid minor
% Interval Density
subplot(3,2,4);
plot([250,650],[250,650],'k','linewidth',2); hold on;errorbar(1000.*CoreIntervalDensity,1000.*RadarIntervalDensity,1000.*sqrt(RadarIntervalDensityVar),'ok','markersize',0.5,'linewidth',2);...
scatter(1000.*CoreIntervalDensity,1000.*RadarIntervalDensity,50,colors,'filled','markeredgecolor',[0.75,0.75,0.75],'linewidth',.5);
title('Interval Density [kgm^{-3}]')
set(gca,'fontsize',12,'fontweight','bold')
axis tight
grid minor
% RMS Average Accumulation
subplot(3,2,5);
plot([0,10],[0,10],'k','linewidth',2); hold on;errorbar(CoreRMSaccumulation,RadarRMSaccumulation,sqrt(RadarRMSaccumulationVar),'ok','markersize',0.5,'linewidth',2);...
scatter(CoreRMSaccumulation,RadarRMSaccumulation,50,colors,'filled','markeredgecolor',[0.75,0.75,0.75],'linewidth',.5);
title('Average Accumulation [m.w.e.]')
set(gca,'fontsize',12,'fontweight','bold')
axis tight
grid minor
% Interval Accumulation
subplot(3,2,6);
plot([0,3],[0,3],'k','linewidth',2); hold on;errorbar(CoreIntervalAccumulation,RadarIntervalAccumulation,sqrt(RadarIntervalAccumulationVar),'ok','markersize',0.5,'linewidth',2);...
scatter(CoreIntervalAccumulation,RadarIntervalAccumulation,50,colors,'filled','markeredgecolor',[0.75,0.75,0.75],'linewidth',.5);
% scatter(CoreIntervalAccumulation,RadarDiffAccumulation,50,colors,'filled','markeredgecolor',[0.75,0.75,0.75],'linewidth',.5);
title('Interval Accumulation [m.w.e.]')
set(gca,'fontsize',12,'fontweight','bold')
axis tight
grid minor
labt = suplabel('Core 7 Validation of Multi-channel Radar Derived Estimates','t');
laby = suplabel('Radar Estimate','y');
labx = suplabel('Core Site 7 Validation','x');
set(labt,'fontsize',20)
set(laby,'fontsize',16,'fontweight','bold')
set(labx,'fontsize',16,'fontweight','bold')