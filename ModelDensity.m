% clear
  addpath '/home/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/Save'
addpath '/home/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/TM'
load('yet_white.mat')

RadarValidIx = find(RadarDepth<=12);
RadDepth = RadarDepth(RadarValidIx);%[0.0697;2.3677;3.3069;4.0175;4.6609;6.5314;8.655;10.1861;11.5721];
RadVrms = RadarRMSv(RadarValidIx);%[0.2357;0.2330;0.2287;0.2249;0.2163;0.2113;0.2122;0.2127;0.2141];
RadDensity = DryCrim(RadVrms);%[322.0;339.3;368.3;394.9;457.7;495.8;490.1;484.9;473.8];
% Impose Herron & Langway Model for Invertable/Stable Vrms Gradient
% Depths for Model Estimate are Solved from Vrms Profile
%         Zhl = [Zlmo;Znmo];
% Zhl = linspace(0.001,30.001,301);Zhl = Zhl(:);
Zhl = linspace(0.001,max(Profile(:,1)),max(Profile(:,1).*10+1));Zhl = Zhl(:);

Hhl = [0.001;diff(Zhl)];
% Assume Mean Annual Temperature at -20C
annualT = -25:1:-15;
% Compute Surface Density from Direct Wave Velocity
surfDensity = 0.3:0.005:0.4;
% Assume Annual Accumulation at 0.5 m.w.e.
annualA = 0.75:-.05:0.1;%0.1:0.05:0.75;
% allocation
% densityHL = cell(length(annualA),length(surfDensity),length(annualT));
% densityAvgHL = cell(length(annualA),length(surfDensity),length(annualT));
% Vhl = cell(length(annualA),length(surfDensity),length(annualT));
densityHL = zeros(length(Zhl),length(annualA)*length(surfDensity)*length(annualT));
densityAvgHL = zeros(length(Zhl),length(annualA)*length(surfDensity)*length(annualT));
Vhl = zeros(length(Zhl),length(annualA)*length(surfDensity)*length(annualT));
iter = 1;
for ii = 1:length(surfDensity)
    for kk = 1:length(annualA);
        for jj = 1:length(annualT)
            
            
            % Compute HL Model Density
            %             densityHL{kk,ii,jj} = herronLangway( Zhl, annualT(jj), surfDensity(ii), annualA(kk));
            densityHL(:,iter) = herronLangway( Zhl, annualT(jj), surfDensity(ii), annualA(kk));
            
            % Compute Cumulative Average Density
            %             cumW = zeros(length(densityHL{kk,ii,jj}),length(densityHL{kk,ii,jj}));
            cumW = zeros(length(densityHL(:,iter)),length(densityHL(:,iter)));
            %             for ll = 1:length(densityHL{kk,ii,jj})
            for ll = 1:length(densityHL(:,iter))
                
                cumW(ll,1:ll) = Hhl(1:ll)./Zhl(ll);
            end
            densityAvgHL(:,iter) = cumW*densityHL(:,iter);
            %             densityAvgHL{kk,ii,jj} = cumW*densityHL{kk,ii,jj};
            %         densityAvgHL = cumsum(densityHL)./(1:length(densityHL))';
            % Estimate Vrms from HL model ~ to Average Density
            Vhl(:,iter) = DryCrimVRMS(densityAvgHL(:,iter));
            %             Vhl{kk,ii,jj} = DryCrimVRMS(densityAvgHL{kk,ii,jj});
            iter = iter+1;
        end
    end
end
        % Find HL Index Nearest the Radar Estimation
        [~,HLix] = min(pdist2(RadDepth,Zhl),[],2);
        
        % Get Average Velocities from HL Estimates at Radar depths
        compZ = Zhl(HLix); compV = Vhl(HLix,:);
        
        % Calculate RMSE of HL and Radar
        RMSE = zeros(1,length(annualA)*length(surfDensity)*length(annualT));
        for ii = 1:length(annualA)*length(surfDensity)*length(annualT)
            RMSE(ii) = sqrt(mean((RadVrms-compV(:,ii)).^2));
        end
        % X = ones(length(Zhl),1)*[1:length(annualA)*length(surfDensity)*length(annualT)];
Z = Zhl*ones(1,length(annualA)*length(surfDensity)*length(annualT));
% X = ones(length(Zhl),1)*[1:length(annualA)*length(surfDensity)*length(annualT)];
X = ones(length(Zhl),1)*densityHL(1,:);
        [~,minIx] = min(RMSE);
        goodIx = find(RMSE < RMSE(minIx)+0.0005 & RMSE > RMSE(minIx)-0.000125);
        QdensityAvgHL = quantile(densityAvgHL(:,goodIx),[.025,.50,0.975],2);
        QdensityHL = quantile(densityHL(:,goodIx),[.025,.50,0.975],2);
        [lowSurfavg] = min(densityAvgHL(1,goodIx));
        [highSurfavg] = max(densityAvgHL(1,goodIx));
        [lowSurf] = min(densityHL(1,goodIx));
        [highSurf] = max(densityHL(1,goodIx));

        BestHLfit = Vhl(:,minIx);
        BestHLdensity = densityHL(:,minIx);
        BestHLavgDensity = DryCrim(BestHLfit);
        MeanHLavgDensity = DryCrim(mean(Vhl(:,goodIx),2));
        varHLavgDensity = var(densityAvgHL(:,goodIx),0,2);
        %%
        % Create To intercept Times from Best Fit Model
        ToHL = (2.*(Zhl))./Vhl(:,minIx);
        RadTo = [mean(LMOto(2,:),2);mean(NMOto,2)];
        % Create Matrix of Radar Velocities
        for jj = 1:length(RadTo)
            Vstack(jj,:) = linspace(RadarRMSv(jj)-2.*sqrt(RadarRMSvar(jj)),RadarRMSv(jj)+2.*sqrt(RadarRMSvar(jj)),500);
        end
        % Velocity Coherence
        Gather = 18900;
        for jj = 1:nChan
            % To Correction
            data(:,jj) = mean(Radar{jj}(1:end,Gather),2);
            directData(:,jj) = mean(directRadar{jj}(1:end,Gather),2);
%             data(:,jj) = mean(Radar{jj}(round(mean([AirTo{Gather}])./dt):end,Gather),2);
%             directData(:,jj) = mean(directRadar{jj}(round(mean([AirTo{Gather}])./dt):end,Gather),2);
        end
        % AGC gain
        scandata = AGCgain(data,250,2);
%         [S,T,V] = GuidedVelocityCoherence(data,ToHL,dt,offsetArray,Vhl(:,goodIx),2,50,1);
%         [S,T,V] = GuidedVelocityCoherence(data,RadTo,dt,offsetArray,RadarRMSv,2,50,1);
%         [S,T,V] = GuidedVelocityCoherence(data,RadTo,dt,offsetArray,Vstack,2,50,1);
        [S,T,V] = VelocityCoherence(scandata(:,:),dt,offsetArray(:,:),.169,.25,1000,2,50,1);
%         clear Tgrid;
        Snorm = S./max(S(:));Tgrid = T(:)*ones(1,1000);
%         Threshold = find(Snorm<0.35);
%         Snorm(Threshold) = 0;
%         Sagc = AGCgain(S,101,2);

%         Snorm = S./max(S(:));Tgrid = T(:)*ones(1,size(Vstack,2));
%         Snorm = S./max(S(:));Tgrid = T(:)*ones(1,length(goodIx));
%         figure();pcolor(V,Tgrid,Snorm);colormap(yet_white);shading interp; axis ij;
        figure();pcolor(V,Tgrid,Snorm);colormap(yet_white);shading interp; axis ij;
        hold on;plot(DirectVelocity{2}(Gather),DirectTo{2}(Gather),'ok');
        plot(DirectVelocity{2}(Gather)-sqrt(DirectVelocityVar{2}(Gather)),DirectTo{2}(Gather),'.k');
        plot(DirectVelocity{2}(Gather)+sqrt(DirectVelocityVar{2}(Gather)),DirectTo{2}(Gather),'.k');
        plot(ReflectionVelocity{1}(Gather),ReflectionTo{1}(Gather),'ok');
        plot(ReflectionVelocity{1}(Gather)-sqrt(ReflectionVelocityVar{1}(Gather)),ReflectionTo{1}(Gather),'.k');
        plot(ReflectionVelocity{1}(Gather)+sqrt(ReflectionVelocityVar{1}(Gather)),ReflectionTo{1}(Gather),'.k');        
        plot(ReflectionVelocity{2}(Gather),ReflectionTo{2}(Gather),'ok')
        plot(ReflectionVelocity{2}(Gather)-sqrt(ReflectionVelocityVar{2}(Gather)),ReflectionTo{2}(Gather),'.k');
        plot(ReflectionVelocity{2}(Gather)+sqrt(ReflectionVelocityVar{2}(Gather)),ReflectionTo{2}(Gather),'.k');        
        plot(ReflectionVelocity{3}(Gather),ReflectionTo{3}(Gather),'ok')
        plot(ReflectionVelocity{3}(Gather)-sqrt(ReflectionVelocityVar{3}(Gather)),ReflectionTo{3}(Gather),'.k');
        plot(ReflectionVelocity{3}(Gather)+sqrt(ReflectionVelocityVar{3}(Gather)),ReflectionTo{3}(Gather),'.k'); 
        plot(ReflectionVelocity{4}(Gather),ReflectionTo{4}(Gather),'ok')
        plot(ReflectionVelocity{4}(Gather)-sqrt(ReflectionVelocityVar{4}(Gather)),ReflectionTo{4}(Gather),'.k');
        plot(ReflectionVelocity{4}(Gather)+sqrt(ReflectionVelocityVar{4}(Gather)),ReflectionTo{4}(Gather),'.k');  
        plot(ReflectionVelocity{5}(Gather),ReflectionTo{5}(Gather),'ok')
        plot(ReflectionVelocity{5}(Gather)-sqrt(ReflectionVelocityVar{5}(Gather)),ReflectionTo{5}(Gather),'.k');
        plot(ReflectionVelocity{5}(Gather)+sqrt(ReflectionVelocityVar{5}(Gather)),ReflectionTo{5}(Gather),'.k');         
        ii = 1;
        HorizonTravelTime = mean(exportReflectionTravelTimes{ii}(:,:,Gather),3);
        DirectTravelTime = mean(exportDirectTravelTimes{ii}(:,:,Gather),3);
        hold on;
        errorbarxycustom(RadarRMSv,RadTo, sqrt(RadarRMSvar), zeros(1,length(RadTo)),'Color','k','LineStyle','None','Linewidth',2);
% herrorbar(RadDensity, RadDepth, sqrt(RadarRMSdensityVar(RadarValidIx)), 'ok')
scatter(RadarRMSv,RadTo,100,colors,'filled','markeredgecolor',[0.75,0.75,0.75],'linewidth',.5);
        [~,sortOffsetIx] = sort(offsetArray);
        dosTiempo = (0:1:size(data(:,1),1)-1)*dt;
        
%         coTo = 124; coVrms = .1912;
%         t8To = 128.5; t8Vrms = .2279;
%         coTo = 113; coVrms = .1908;
%         t8To = 108; t8Vrms = .214;
%         coTo = 91.6; coVrms = .1905;
%         coWindow = round(coTo./dt)-15:round(coTo./dt)+15;
%         t8To = 95.8; t8Vrms = .2124;
%         t8Window = round(t8To./dt)-15:round(t8To./dt)+15;
%         coTo = 42.8; coVrms = .2255;
%         coWindow = round(coTo./dt)-15:round(coTo./dt)+15;
%         t8To = 42.8; t8Vrms = .2149;
%         t8Window = round(t8To./dt)-15:round(t8To./dt)+15;
        coTo = 38.8; coVrms = .233;
        coWindow = round(coTo./dt)-15:round(coTo./dt)+15;
        t8To = 35.4; t8Vrms = .2235;
        t8Window = round(t8To./dt)-15:round(t8To./dt)+15;

        % NMO Correction
        [coNMO] = NMO(data,dt,offsetArray,coTo,coVrms,100);
        [t8NMO] = NMO(data,dt,offsetArray,t8To,t8Vrms,100);
        
        dude = NMO(directData,dt,offsetArray,0,0.299,100);
        
%         ToIntHL = diff(toHL(HLix));
% [~,minsortIx] = sort(densityHL(1,:));
% sortDensityHL = densityHL(:,minsortIx);
% sortDensityHL = densityHL(:,sortmeanDensityIx);
%% Wiggle Plot if isWiggle
% for ii = 11%1:selection
% Open CMP Arrival Pick Wiggle
        % Open Figure
%         figstr = ['Traversed CMP Reflection Gather' num2str(ii,'%02d')];
%         h = figure('Name',figstr,'NumberTitle','off');
%        set(h, 'Position', [500 500 600 850])
        % Create Wiggle Plot
        figure();
        wiggle(dosTiempo,offsetArray,data)
% imagesc(offsetArray(sortOffsetIx),dosTiempo,data(:,sortOffsetIx));colormap(LateNite)
        hold on
        for jj = 1:size(HorizonTravelTime,2)
%         plot(offsetArray(sortOffsetIx),HorizonTravelTime(sortOffsetIx,jj)+AirTo{Gather}.*ones(length(offsetArray),1),'r','linewidth',1)
%         plot(offsetArray(sortOffsetIx),HorizonTravelTime(sortOffsetIx,jj)+AirTo{Gather}.*ones(length(offsetArray),1),'.k','linewidth',1)
        plot(offsetArray(sortOffsetIx),HorizonTravelTime(sortOffsetIx,jj),'r','linewidth',1)
        plot(offsetArray(sortOffsetIx),HorizonTravelTime(sortOffsetIx,jj),'xk','linewidth',2)
        end
%         plot(offsetArray(sortOffsetIx),synTime+AirTo{Gather}.*ones(1,length(offsetArray)),'b')
%         set(gca,'LineWidth',2,'FontSize',14, 'FontWeight'...
%         ,'bold')
% %         titlestr = ['Radar Reflection Gather Location ' num2str(CMP(selection(ii)).Location,'%5.2f')];
% %         title(titlestr,'FontWeight','Bold');
%         set(gca,'XTick',[0 4 8 12], 'XTickLabel',[0 4 8 12])
%         xlabel('Offset [m]')
%         ylabel('Two-Way Time [ns]')
%         
%         figure();
%         wiggle(dosTiempo,offsetArray,data)
% imagesc(offsetArray(sortOffsetIx),dosTiempo,data(:,sortOffsetIx));colormap(LateNite)
        hold on
        plot(offsetArray(sortOffsetIx),offsetArray(sortOffsetIx)./0.3,'b')
        for jj = 1:size(DirectTravelTime,2)
        plot(offsetArray(sortOffsetIx),DirectTravelTime(sortOffsetIx,jj),'r','linewidth',1)
        plot(offsetArray(sortOffsetIx),DirectTravelTime(sortOffsetIx,jj),'xk','linewidth',2)
        end
%         wiggle(t,CMP(ii).Offset,CMP(ii).Data)
%         wiggle(time{1},CMP(selection(ii)).Offset,filtTimeCMP{ii})
        set(gca,'LineWidth',2,'FontSize',14, 'FontWeight'...
        ,'bold')
%         titlestr = ['Radar Reflection Gather Location ' num2str(CMP(selection(ii)).Location,'%5.2f')];
%         title(titlestr,'FontWeight','Bold');
        set(gca,'XTick',[0 4 8 12], 'XTickLabel',[0 4 8 12])
        xlabel('Offset [m]')
        ylabel('Two-Way Time [ns]')
        
%         figure();
%         subplot(2,1,1)
%         wiggle(dosTiempo(coWindow),offsetArray,0.025.*coNMO(coWindow,:))
%         title(['To = ',num2str(coTo),' V = ',num2str(coVrms),' from Coherence Scan'])
%         subplot(2,1,2)
%         wiggle(dosTiempo(t8Window),offsetArray,0.025.*t8NMO(t8Window,:))
%         title(['To = ',num2str(t8To),' V = ',num2str(t8Vrms),' from t^{2}-x^{2} Regression'])
        
% end


%% Plots
Heavy = 0;
if Heavy
        solX = X(1,minIx);
        q025X = QdensityHL(1,1);
        q975X = QdensityHL(1,3);
        
figure();
h1 = subplot(2,2,1);
pcolor(X,Z,densityAvgHL);shading flat;axis ij;colormap(flipud(bone));freezeColors;hold on; 
xtick = [surfDensity(2),(linspace(lowSurfavg,highSurfavg,3)),surfDensity(end)];
set(h1,'XAxisLocation','Top','XTick',xtick);
pcolor(X(:,goodIx),Z(:,goodIx),densityAvgHL(:,goodIx));colormap(yet_white);shading flat;axis ij;
h = colorbar;
ylabel(h,'Density [kgm^{-3}]','rotation',270,'fontsize',14)
hy = get(h,'YLabel');
set(hy, 'Units','Normalized')
pos = get(hy,'Position');
set(hy,'Position',pos.*[1.5,1,1])
caxis([min(min(densityAvgHL(:,goodIx))),max(max(densityAvgHL(:,goodIx)))])
freezeColors;
plot(ones(length(Zhl),1).*solX,Zhl,'k','linewidth',2)
plot(ones(length(Zhl),1).*q025X,Zhl,'--k','linewidth',1)
plot(ones(length(Zhl),1).*q975X,Zhl,'--k','linewidth',1)
% pcolor(X(:,solIx),Z(:,solIx),densityAvgHL(:,solIx));colormap(black);shading flat;axis ij;
% freezeColors;

set( h, 'YDir', 'reverse' );
% set(h1,'XTickLabel',{'0.3','0.35','0.365','0.38','0.4'})
h1.XAxis.TickLabelFormat = '%.3f';
set(gca,'fontsize',12,'fontweight','bold')
ylabel('Depth [m]')
xlabel('Herron-Langway Average Density','fontsize',16)
% colorbar
subplot(2,2,2)
hold on;
% stairs(QdensityAvgHL(:,2),Zhl,'b','linewidth',1);
% stairs(QdensityAvgHL(:,1),Zhl,'--k','linewidth',1);
% stairs(QdensityAvgHL(:,3),Zhl,'--k','linewidth',1);
stairs(BestHLavgDensity,Zhl,'k','linewidth',2);hold on; 
% stairs(densityAvgHL(:,goodIx(lowEndIx)),Zhl,'-- k','linewidth',1);
% stairs(densityAvgHL(:,goodIx(highEndIx)),Zhl,'-- k','linewidth',1);
% stairs(MeanHLavgDensity,Zhl,'b','linewidth',2);
% stairs(MeanHLavgDensity-sqrt(varHLavgDensity),Zhl,'--b','linewidth',1);
% stairs(MeanHLavgDensity+sqrt(varHLavgDensity),Zhl,'--b','linewidth',1);
% stairs(MeanHLavgDensity-2.*sqrt(varHLavgDensity),Zhl,'--b','linewidth',1);
% stairs(MeanHLavgDensity+2.*sqrt(varHLavgDensity),Zhl,'--b','linewidth',1);
stairs(Profile(:,4),Profile(:,1),'r','linewidth',2);axis ij;axis tight
% plot(RadDensity,RadDepth,'ok','linewidth',2,'markersize',10)
% errorbarxycustom(RadDensity, RadDepth, sqrt(RadarRMSdensityVar(RadarValidIx)), sqrt(RadarDepthVar(RadarValidIx)),'Color','k','LineStyle','None','Linewidth',2);
% herrorbar(RadDensity, RadDepth, sqrt(RadarRMSdensityVar(RadarValidIx)), 'ok')
scatter(RadDensity,RadDepth,100,colors(RadarValidIx,:),'filled','markeredgecolor',[0.75,0.75,0.75],'linewidth',.5);
legend('Herron-Langway Fit','Core 7 Log','Radar Estimate','location','southwest')
stairs(QdensityAvgHL(:,1),Zhl,'--k','linewidth',1);
stairs(QdensityAvgHL(:,3),Zhl,'--k','linewidth',1);
stairs(BestHLavgDensity,Zhl,'k','linewidth',2);hold on; 
stairs(Profile(:,4),Profile(:,1),'r','linewidth',2);axis ij;axis tight
errorbarxycustom(RadDensity, RadDepth, sqrt(RadarRMSdensityVar(RadarValidIx)), sqrt(RadarDepthVar(RadarValidIx)),'Color','k','LineStyle','None','Linewidth',2);
% herrorbar(RadDensity, RadDepth, sqrt(RadarRMSdensityVar(RadarValidIx)), 'ok')
scatter(RadDensity,RadDepth,100,colors(RadarValidIx,:),'filled','markeredgecolor',[0.75,0.75,0.75],'linewidth',.5);
set(gca,'fontsize',12,'fontweight','bold')
xlabel('Density [kgm^{-3}]')
title('Radar Informed Solutions','fontsize',16)
% figure();
h3 = subplot(2,2,3);
pcolor(X,Z,densityHL);shading flat;axis ij;colormap(flipud(bone));freezeColors;hold on;
xtick = [surfDensity(2),(linspace(lowSurf,highSurf,3)),surfDensity(end)];
set(h3,'XAxisLocation','Top','XTick',xtick);
pcolor(X(:,goodIx),Z(:,goodIx),densityHL(:,goodIx));colormap(yet_white);shading flat;axis ij;
h = colorbar;
ylabel(h,'Density [kgm^{-3}]','rotation',270,'fontsize',14)
hy = get(h,'YLabel');
set(hy, 'Units','Normalized')
pos = get(hy,'Position');
set(hy,'Position',pos.*[1.75,1,1])
caxis([min(min(densityHL(:,goodIx))),max(max(densityHL(:,goodIx)))]);
freezeColors;
plot(ones(length(Zhl),1).*solX,Zhl,'k','linewidth',2)
plot(ones(length(Zhl),1).*q025X,Zhl,'--k','linewidth',1)
plot(ones(length(Zhl),1).*q975X,Zhl,'--k','linewidth',1)
% pcolor(X(:,solIx),Z(:,solIx),densityHL(:,solIx));colormap(black);shading flat;axis ij;
% freezeColors;

set( h, 'YDir', 'reverse' );
% set(h3,'XTickLabel',{'0.3','0.35','0.365','0.38','0.4'})
h3.XAxis.TickLabelFormat = '%.3f';

set(gca,'fontsize',12,'fontweight','bold')
ylabel('Depth [m]')
xlabel('Herron-Langway Density','fontsize',16)
% colorbar
subplot(2,2,4)
hold on;
% stairs(QdensityHL(:,2),Zhl,'b','linewidth',1);

stairs(BestHLdensity,Zhl,'k','linewidth',2);hold on; 
% stairs(densityHL(:,goodIx(lowEndIx)),Zhl,'-- k','linewidth',1);
% stairs(densityHL(:,goodIx(highEndIx)),Zhl,'-- k','linewidth',1);
stairs(Profile(:,3),Profile(:,1),'r','linewidth',2);axis ij;axis tight
% errorbarxycustom(RadarIntervalDensity(RadarValidIx), RadDepth, sqrt(RadarIntervalDensityVar(RadarValidIx)), sqrt(RadarDepthVar(RadarValidIx)),'Color','k','LineStyle','none','Linewidth',2);
% herrorbar(RadarIntervalDensity(RadarValidIx), RadDepth, sqrt(RadarIntervalDensityVar(RadarValidIx)),'ok')
scatter(RadarIntervalDensity(RadarValidIx),RadDepth,100,colors(RadarValidIx,:),'filled','markeredgecolor',[0.75,0.75,0.75],'linewidth',.5);
legend('Herron-Langway Fit','Core 7 Log','Radar Estimate','location','southwest')
% cla
stairs(QdensityHL(:,1),Zhl,'--k','linewidth',1);
stairs(QdensityHL(:,3),Zhl,'--k','linewidth',1);
stairs(BestHLdensity,Zhl,'k','linewidth',2);hold on; 
% stairs(densityHL(:,goodIx(lowEndIx)),Zhl,'-- k','linewidth',1);
% stairs(densityHL(:,goodIx(highEndIx)),Zhl,'-- k','linewidth',1);
stairs(Profile(:,3),Profile(:,1),'r','linewidth',2);axis ij;axis tight
errorbarxycustom(RadarIntervalDensity(RadarValidIx), RadDepth, sqrt(RadarIntervalDensityVar(RadarValidIx)), sqrt(RadarDepthVar(RadarValidIx)),'Color','k','LineStyle','none','Linewidth',2);
% herrorbar(RadarIntervalDensity(RadarValidIx), RadDepth, sqrt(RadarIntervalDensityVar(RadarValidIx)),'ok')
scatter(RadarIntervalDensity(RadarValidIx),RadDepth,100,colors(RadarValidIx,:),'filled','markeredgecolor',[0.75,0.75,0.75],'linewidth',.5);
set(gca,'fontsize',12,'fontweight','bold')
xlabel('Density [kgm^{-3}]')
title('Radar Informed Solutions','fontsize',16)
%%
% clear
end