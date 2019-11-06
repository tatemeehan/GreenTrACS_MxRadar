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
        figure();imagesc(Traverse{ii}./1000,DepthAxis{ii},1000.*DensityModel{ii});colormap(yet_white);freezeColors;hold on;
        imagesc(Traverse{ii}./1000,DepthAxis{ii},sign(RadarDepth{ii}),'AlphaData',WiggleAlpha);colormap([0,0,0]);freezeColors;
        colormap(yet_white);hlay = colorbar; set(hlay,'YDir','reverse','fontsize',14,'fontweight','bold','Ticks',[310,350,400,450,500,550,600,625]);
        set(get(hlay,'ylabel'),'String','Density [kg/m^{3}]', 'rotation', 270,'Units', 'Normalized', 'Position', [4, 0.5, 0])
        title('Density Tomogram')
        xlabel('Distance (km)')
        ylabel('Depth (m)','rotation',270, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0])
        set(gca,'fontsize',14,'fontweight','bold')
        set(gca,'YTick',[0,2.5,5,7.5,10,12.5,15,17.5,20,22.5])
        
        % Plot Density Anamoly, Overlay Peak Amplitudes
%         figure();imagesc(Traverse{ii}./1000,DepthAxis{ii},1000.*DensityAnomalyModel{ii});colormap(colorbrew);caxis([-15,15]);freezeColors;hold on;
% %         plot(Traverse{ii},depth550,'--k','linewidth',3);
%         imagesc(Traverse{ii},DepthAxis{ii},sign(RadarDepth{ii}),'AlphaData',WiggleAlpha);colormap([0,0,0]);freezeColors;
% %         colormap(SplitJet);hlay = colorbar; %set(hlay,'YDir','reverse','fontsize',14,'fontweight','bold');
%         colormap(colorbrew);hlay = colorbar;
        figure();
        imagesc(Traverse{ii}./1000,DepthAxis{ii},1000.*DensityAnomalyModel{ii});colormap(colorbrew);caxis([-15,15]);freezeColors;hold on;
        imagesc(Traverse{ii}./1000,DepthAxis{ii},sign(RadarDepth{ii}),'AlphaData',WiggleAlpha);colormap([0,0,0]);freezeColors;
        colormap(colorbrew);hlay = colorbar;
        set(get(hlay,'ylabel'),'String','Deviation from Mean Density [kg/m^{3}]', 'rotation', 270,'Units', 'Normalized', 'Position', [4, 0.5, 0])
        title('Density Anomaly')
        xlabel('Distance (km)')
        ylabel('Depth (m)','rotation',270, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0])
        set(gca,'fontsize',14,'fontweight','bold')
        set(gca,'YTick',[0,2.5,5,7.5,10,12.5,15,17.5,20,22.5])
        
        figure();
        hDA1 = subplot(2,1,1);
        % Plot Mean Average Density Deviation
%         plot(Traverse{ii},1000.*MeanDensityDeviation{ii},'k','linewidth',3)
%         title('Mean Average Deviation - Density [kg/m^{3}]')
        % Plot Surface Density Deviation
        plot(Traverse{ii}./1000,1000.*SurfaceDensityDeviation{ii},'k','linewidth',3)
        title('Surface Density - Deviation from Mean (kg/m^{3})')
        grid on
        set(gca,'fontsize',14,'fontweight','bold')
%         set(hDA1,'units','normalized')
%         hDA1pos = get(hDA1,'Position');
%         set(hDA1,'position',[hDA1pos(1),hDA1pos(2)+.25, hDA1pos(3), hDA1pos(2)-.25])
        
        hDA2 = subplot(2,1,2);
        % Plot Density Anamoly, Overlay Peak Amplitudes
        imagesc(Traverse{ii}./1000,DepthAxis{ii},1000.*DensityAnomalyModel{ii});colormap(colorbrew);caxis([-15,15]);freezeColors;hold on;
        imagesc(Traverse{ii}./1000,DepthAxis{ii},sign(RadarDepth{ii}),'AlphaData',WiggleAlpha);colormap([0,0,0]);freezeColors;
        colormap(colorbrew);hlay = colorbar;
%         colormap(SplitJet);hlay = colorbar; %set(hlay,'YDir','reverse','fontsize',14,'fontweight','bold');
        set(get(hlay,'ylabel'),'String','Deviation from Mean Density (kg/m^{3})', 'rotation', 270,'Units', 'Normalized', 'Position', [4, 0.5, 0])
        title('Density Anomaly')
        xlabel('Distance (km)')
        ylabel('Depth (m)','rotation',270, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0])
        set(gca,'fontsize',14,'fontweight','bold')
        set(gca,'YTick',[0,2.5,5,7.5,10,12.5,15,17.5,20,22.5])
        
        % Plot Depth-Age Tomography
        figure();imagesc(Traverse{ii}./1000,DepthAxis{ii},AgeModel{ii});colormap(yet_white);freezeColors;hold on;
        imagesc(Traverse{ii}./1000,DepthAxis{ii},sign(RadarDepth{ii}),'AlphaData',WiggleAlpha);colormap([0,0,0]);freezeColors;
        colormap(yet_white);hlay = colorbar; set(hlay,'YDir','reverse','fontsize',14,'fontweight','bold','Ticks',[Ages{ii}]);
        for kk = 1:size(DepthAge{ii},2)
        plot(linspace(Traverse{ii}(1)./1000,Traverse{ii}(end)./1000,length(Traverse{ii})),DepthAge{ii}(:,kk),'k','linewidth',3)
        end
        set(get(hlay,'ylabel'),'String','Age (a)', 'rotation', 270,'Units', 'Normalized', 'Position', [4, 0.5, 0])
        set(gca,'YTick',[0,2.5,5,7.5,10,12.5,15,17.5,20,22.5])
        title('Isochronogram')
        xlabel('Distance (km)')
        ylabel('Depth (m)','rotation',270, 'Units', 'Normalized', 'Position', [-0.08, 0.5, 0])
        set(gca,'fontsize',14,'fontweight','bold')            
  
    end
    %% Image Depth Section
    for ii = 1:nFiles
%         figure();imagesc(Traverse{ii},DepthAxis{ii},RadarDepth{ii});
%         figure();imagesc(Traverse{ii}./1000,DepthAxis{ii},AGCgain(RadarDepth{ii},size(RadarDepth{ii},1)./round(5),2));
%         colormap(cmapAdapt(RadarDepth{ii},colorbrew));hold on;
 figure();imagesc(Traverse{ii}./1000,DepositionAxis{ii},AGCgain(RadarDeposition{ii},size(RadarDeposition{ii},1)./round(5),2));
        colormap(cmapAdapt(RadarDeposition{ii},colorbrew));hold on;
%         for kk = 1:size(DepthAge{ii},2)
%             plot(linspace(Traverse{ii}(1)./1000,Traverse{ii}(end)./1000,length(Traverse{ii})),DepthAge{ii}(:,kk),'color',c3,'linewidth',2)
%         end
%         for kk = 1:size(DepthAge{ii},2)
%         plot(Traverse{ii},DepthAge{ii}(:,kk),'k','linewidth',3)
%         end
%         title('Core 15 Spur West - Depth Section')
        title('Core 15 Spur West - Age Section')

        xlabel('Distance (km)')
%         ylabel('Depth (m)','rotation',270, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0])
        ylabel('age (a)','rotation',270, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0])
        set(gca,'fontsize',14,'fontweight','bold')
        

        
        % Plot RadarGram with Isochrones
%         figure();imagesc(Traverse{ii},DepthAxis{ii},RadarDepth{ii});
        figure();imagesc(Traverse{ii}./1000,DepthAxis{ii},AGCgain(RadarDepth{ii},size(RadarDepth{ii},1)./round(3.5),2));
        colormap(cmapAdapt(RadarDepth{ii},colorbrew));hold on;
        for kk = 1:size(DepthAge{ii},2)
        plot(linspace(Traverse{ii}(1)./1000,Traverse{ii}(end)./1000,length(Traverse{ii})),DepthAge{ii}(:,kk),'color',c3,'linewidth',2)
        end
        
        % Compare Time to Depth Images
        compareIx = 7550:11500;
        compareIy = 165:565;
        % Salt and Pepper Time Image
        figure();
             pcolor(Traverse{ii}(compareIx)./1000,tStack{ii}(compareIy,1),AGCgain(RadarStack{1}(compareIy,compareIx),size(Radar{ii}(:,:),1)./round(7.5),2)); shading interp; axis ij
        colormap(cmapAdapt(RadarStack{1}(compareIy,compareIx),colorbrew));hold on;
%         subplot(2,1,1)
%      imagesc(Traverse{ii}(:)./1000,tStack{ii}(:,1),AGCgain(RadarNMO{8}(:,:),size(Radar{ii}(:,:),1)./round(3.5),2));
%         colormap(cmapAdapt(Radar{4}(compareIy,compareIx),colorbrew));hold on;
%         imagesc(Traverse{ii}(compareIx)./1000,tStack{ii}(compareIy,1),AGCgain(RadarNMO{8}(compareIy,compareIx),size(Radar{ii}(compareIy,compareIx),1)./round(3.5),2));
%         colormap(cmapAdapt(Radar{4}(compareIy,compareIx),colorbrew));hold on;
%         title('Core 15 Spur West - Time Section')
        title('Time Section')
        xlabel('Distance (km)')
        ylabel('Travel-Time (ns)','rotation',270, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0])
        set(gca,'fontsize',14,'fontweight','bold','layer','top')
                daspect([1 9.5 1])

%         axis square
%         compareIy = 150:550;
                compareIy = 150:1000;

        % Bread and Butter Depth Image
        figure();
%         subplot(2,1,2)
%         imagesc(Traverse{ii}(compareIx)./1000,DepositionAxis{ii}(compareIy),AGCgain(RadarDeposition{ii}(compareIy,compareIx),size(RadarDeposition{ii}(compareIy,compareIx),1)./round(3.5),2));
        pcolor(Traverse{ii}(compareIx)./1000,DepositionAxis{ii}(compareIy),AGCgain(RadarDeposition{ii}(compareIy,compareIx),size(RadarDeposition{ii}(compareIy,compareIx),1)./round(5),2));shading interp, axis ij
        colormap(cmapAdapt(RadarDeposition{ii}(compareIy,compareIx),colorbrew));hold on;
                
%                 for kk = 1:length(isochronePick)
%                     plot(Traverse{ii}(compareIx)./1000,isochronePick{kk}(compareIx),'k','linewidth',2)
%                 end

%         imagesc(Traverse{ii}(compareIx)./1000,DepthAxis{ii}(compareIy),AGCgain(RadarDepth{ii}(compareIy,compareIx),size(RadarDepth{ii}(compareIy,compareIx),1)./round(3.5),2));
%         colormap(cmapAdapt(RadarDepth{ii}(compareIy,compareIx),colorbrew));hold on;
%         title('Core 15 Spur West - Depth Section')
%         title('Core 15 Spur West - Age Section')

%         title('Depth Section')
        title('Age Section')
        xlabel('Distance (km)')
%         ylabel('Depth (m)','rotation',270, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0])
        ylabel('Age (a)','rotation',270, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0])

        set(gca,'fontsize',14,'fontweight','bold')
%         daspect([1 1 1])
        daspect([1 2 1])

%         text(.65,-.25,'1000x Vertical Exaggeration','units','normalized','fontsize',10,'fontweight','bold')
%         axis square

        figure();
        compareIy = 150:600;
%         subplot(2,1,2)
%         imagesc(Traverse{ii}(compareIx)./1000,DepthAxis{ii}(compareIy),AGCgain(RadarDepth{ii}(compareIy,compareIx),size(RadarDepth{ii}(compareIy,compareIx),1)./round(3.5),2));
%         colormap(cmapAdapt(RadarDepth{ii}(compareIy,compareIx),colorbrew));hold on;
        pcolor(Traverse{ii}(compareIx)./1000,DepthAxis{ii}(compareIy),AGCgain(RadarDepth{ii}(compareIy,compareIx),size(RadarDepth{ii}(compareIy,compareIx),1)./round(3.5),2)); shading interp; axis ij;
        colormap(cmapAdapt(RadarDepth{ii}(compareIy,compareIx),colorbrew));hold on;
%         title('Core 15 Spur West - Depth Section')
        title('Depth Section')
        xlabel('Distance (km)')
        ylabel('Depth (m)','rotation',270, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0])

        set(gca,'fontsize',14,'fontweight','bold','layer','top')
        daspect([1 1 1])

        text(.65,-.25,'1000x Vertical Exaggeration','units','normalized','fontsize',10,'fontweight','bold')
%         axis square
 
        % Plot Depth-Age Tomography
        figure();pcolor(Traverse{ii}(compareIx)./1000,DepthAxis{ii}(compareIy),bestAgeModel{ii}(compareIy,compareIx));colormap(yet_black);freezeColors;hold on;
%         imagesc(Traverse{ii},DepthAxis{ii},sign(RadarDepth{ii}),'AlphaData',WiggleAlpha);colormap([0,0,0]);freezeColors;
        colormap(yet_white); shading interp;axis ij
        %hlay = colorbar; set(hlay,'YDir','reverse','fontsize',14,'fontweight','bold','Ticks',[5.5:5:20.5],'TickLabels',[2012,2007,2002,1997]);
%         for kk = 1:size(DepthAge{ii},2)
%         plot(Traverse{ii}(compareIx)./1000,DepthAge{ii}(compareIx,kk),'k','linewidth',1.5)
%         end
          contour(Traverse{ii}(compareIx)./1000,DepthAxis{ii}(compareIy),bestAgeModel{ii}(compareIy,compareIx),7:5:17,'k','linewidth',1.5)
%         set(get(hlay,'ylabel'),'String','Age (a)', 'rotation', 270,'Units', 'Normalized', 'Position', [4, 0.5, 0])
%         set(gca,'YTick',[4,6,8,10])
%         title('Core 15 Spur West - Age-Depth Section')
                title('Age-Depth Section')
        xlabel('Distance (km)')
        ylabel('Depth (m)','rotation',270, 'Units', 'Normalized', 'Position', [-0.08, 0.5, 0])
        set(gca,'fontsize',14,'fontweight','bold','layer','top') 
        daspect([1 1 1])
        text(.65,-.25,'1000x Vertical Exaggeration','units','normalized','fontsize',10,'fontweight','bold')
        
        % Age-Travel-Time Model
        compareIx = 7550:11500;
        compareIy = 165:565;
        figure();pcolor(Traverse{ii}(compareIx)./1000,TimeAxis{ii}(compareIy),pseudoAgeModel{ii}(compareIy,compareIx));
        hold on;shading interp;axis ij;
        colormap(yet_white);
        tmpAges = 7:5:17;
        for kk = 1:length(tmpAges)
        [~,tmpIx] = min(abs(pseudoAgeModel{ii}(compareIy,compareIx)-tmpAges(kk)));
        psuedoAgeIx(kk,:) = tmpIx;
        end
%         psuedoAgeIx = find(pseudo
        contour(Traverse{ii}(compareIx)./1000,TimeAxis{ii}(compareIy),pseudoAgeModel{ii}(compareIy,compareIx),7:5:17,'k','linewidth',1.5)
 
        title('Age-Travel-Time Section')
        xlabel('Distance (km)')
        ylabel('Travel-Time (ns)','rotation',270, 'Units', 'Normalized', 'Position', [-0.08, 0.5, 0])
        set(gca,'fontsize',14,'fontweight','bold','layer','top') 
        daspect([1 8.5 1])
        
        % Plot Age-Travel-Time Perturbations
        tmpResidual = updatePseudoAgeModel{ii}(compareIy,compareIx)-pseudoAgeModel{ii}(compareIy,compareIx);
        figure();
        pcolor(Traverse{ii}(compareIx)./1000,TimeAxis{ii}(compareIy),tmpResidual);
        shading interp; axis ij; colormap(cmapAdapt(tmpResidual,SplitJet));
        freezeColors;
        c = colorbar; c.Location = 'northoutside';c.Label.String = 'Age Perturbations (\Delta a)';
        c.FontSize = 12; c.Label.FontSize = 16;
         set(gca,'fontsize',14,'fontweight','bold','layer','top') 
        daspect([1 8.5 1])
                xlabel('Distance (km)')
        ylabel('Travel-Time (ns)','rotation',270, 'Units', 'Normalized', 'Position', [-0.08, 0.5, 0])
        
                % Plot Recalculated Age-Travel-Time Perturbations
                compareIy = 150:600;
        tmpResidual = perturbations(compareIy,compareIx);%updatePseudoAgeModel{ii}(compareIy,compareIx)-pseudoAgeModel{ii}(compareIy,compareIx);
        figure();
        pcolor(Traverse{ii}(compareIx)./1000,DepthAxis{ii}(compareIy),tmpResidual);
        shading interp; axis ij;%caxis([-0.5,0.5])% colormap(cmapAdapt(tmpResidual,SplitJet));
        colormap(cmapAdapt(tmpResidual,SplitJet))
%         freezeColors;
        c = colorbar; c.Location = 'northoutside';c.Label.String = 'Age Perturbations (\Delta a)';
        c.FontSize = 12; c.Label.FontSize = 16;
         set(gca,'fontsize',14,'fontweight','bold','layer','top') 
        daspect([1 1 1])
                xlabel('Distance (km)')
        ylabel('Depth (m)','rotation',270, 'Units', 'Normalized', 'Position', [-0.08, 0.5, 0])
                text(.65,-.25,'1000x Vertical Exaggeration','units','normalized','fontsize',10,'fontweight','bold')

    end
    
    % Plot Accumulation
    q = quantile(tmpAccum(:),[.025,0.975]);
    figure();
%     subplot(4,1,1)
    plot(Traverse{ii}./1000,AverageAccumulation{ii},'k','linewidth',2);
    hold on;
    plot(Traverse{ii}./1000,AverageAccumulation2{ii},'r','linewidth',2)
    plot(Traverse{ii}./1000,AverageAccumulation3{ii},'b','linewidth',2)
    set(gca,'fontsize',10,'fontweight','bold')

    figure();
    subplot(3,1,1)
    % Initial Accumulation
    imagesc(Traverse{ii}./1000,DepthAxis{ii},initAccum,q);colormap(yet_white);
    set(gca,'fontsize',10,'fontweight','bold','xticklabel',[])
    
    subplot(3,1,2)
    % First Update Accumulation
    pcolor(Traverse{ii}./1000,DepthAxis{ii},instantDoom);shading interp;colormap(yet_white); axis ij;
    set(gca,'fontsize',10,'fontweight','bold','layer','top','xticklabel',[]);caxis([q]);

    subplot(3,1,3)
    % Sencond Update Accumulation
    pcolor(Traverse{ii}./1000,DepthAxis{ii},tmpAccum);shading interp;colormap(yet_white);axis ij;
    set(gca,'fontsize',10,'fontweight','bold','layer','top');caxis([q])

    % Compare GTC15 SMB and MXHL SMB
    x = Annuals;
X = [x;x];%flipud(x)]
X = sort(X(:));
Xx = [X(2:end);flipud(X(2:end-2))];
y = bin(:,1);y = y(:);
% Stairs
X = X(2:end-1);
Y = [y,y]'; Y = Y(:);
ci = sqrt(binVar(:,1));
yl = y - ci;
yh = y+ci;
% stairs(x,yl);stairs(x,yh)
Yyh = [yh,yh]';Yyh = Yyh(:);
Yyl = [flipud(yl),flipud(yl)]';Yyl = Yyl(:);
Yy = [Yyh;Yyl];
    figure();
    stairs([Accumulation(1,15);Accumulation(1:33,15)],[2017;Accumulation(1:33,17)],'k','linewidth',2.5);hold on
% patch(Yy,Xx,[0.75 0 0],'FaceAlpha',.25); hold on
plot(Y,X,'color',[0.75,0,0],'linewidth',2.5)
ylabel('Year (CE)')
xlabel('SMB (m w.e. a^{-1})')
title({['Annual Surface Mass Balance']; ['Jan.1984 - Jan.2017']})
set(gca,'fontweight','bold','fontsize',14)
legend('GTC15','MXHL','location','southeast')
ylim([1984 2017])
yticks(fliplr([2017,2010,2000,1990,1984]))

figure();
stairs(GTCdepthAge(:,2),GTCdepthAge(:,1),'k')
hold on;
tmpAgeModel = abs(bestAgeModel{1}(:,1)-(str2num(Year{1})+dayofyear/365));
stairs(tmpAgeModel,DepthAxis{1},'r');
% %     Zstair = [DepthAxis{1};DepthAxis{1}];
% Zstair = [Annuals(1:end-1),flipud(Annuals(1:end-1))];
%     [Zstair,Zix] = sort(Zstair);
%     Zstairz = [Zstair(1:end);flipud(Zstair(1:end))]';
%     Zstairz = flipud(Zstairz);
%     mxhlSMBlow = [bin(:,1)-sqrt(binVar(:,1));bin(:,1)-sqrt(binVar(:,1))];
%     mxhlSMBlow = mxhlSMBlow(Zix);
%     mxhlSMBhi = [bin(:,1)+sqrt(binVar(:,1));bin(:,1)+sqrt(binVar(:,1))];
%     mxhlSMBhi = mxhlSMBhi(Zix);
%     mxhlSMBci = [mxhlSMBlow(:,:),flipud(mxhlSMBhi(:,:))];
%     
%     figure();
%     patch(mxhlSMBci(2:2:end,1),Zstairz(1:2:end-1,1),[.5 0 0 ],'FaceAlpha',0.35);hold on;
%     stairs(bin(:,1),Annuals(1:end-1));
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
            shadedErrorBarT8(Traverse{ii}./1000,SnowWaterEqv{rh,ii},...
                sqrt(SnowWaterEqvVar{rh,ii}),1,{'Color',c1,'linewidth',1.5});
            hold on;
        end
        freezeColors
        axis tight
        axis ij
        grid on
        title('Average Annual Accumulation [m w.e. a^{-1}]')
        subplot(3,1,2)
        for rh = 1:nReflectionHorizon
            shadedErrorBarT8(Traverse{ii}./1000,Density{rh,ii}.*1000,...
                sqrt(DensityVar{rh,ii}).*1000,1,{'Color',[0,0,0],'linewidth',1.5});
            hold on;
        end
        freezeColors
        axis tight
        axis ij
        grid on
        title('Average Snow Density [kg/m^{3}]')
        subplot(3,1,3)
        for rh = 1:nReflectionHorizon
            shadedErrorBarT8(Traverse{ii}./1000,Depth{rh,ii},...
                sqrt(DepthVar{rh,ii}),1,{'Color',c3,'linewidth',1.5});
            hold on;
        end
        freezeColors
        axis tight
        axis ij
        grid on
        title('Snow Depth [m]')
        xlabel('Distance [km]')
        set(findobj(gcf,'type','axes'),'FontName','FreeSerif','FontSize',12,...
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
%             shadedErrorBarT8(distance,dhDepth{dh,ii},...
%                 sqrt(dhDepthVar{dh,ii}),1,{'Color',[1,0.81,0],'linewidth',1.5});
%             hold on;
%         end
        for rh = 1:nReflectionHorizon
            shadedErrorBarT8(distance./1000,Depth{rh,ii},...
                sqrt(DepthVar{rh,ii}),1,{'Color',[0,0,0],'linewidth',1.5});

%                 sqrt(DepthVar{rh,ii}),1,{'Color',[1,0.81,0],'linewidth',1.5});

            hold on;
        end
        freezeColors
        axis ij
        axis tight
        grid on
        ylim([1.6 2.2])
        set(gca,'ytick',[1.6,1.8,2.0,2.2])
        set(gca,'xticklabel',[])
        title('Snow Depth (m)')
        subplot(3,1,2)
        for dh = 2:nDirectHorizon
            shadedErrorBarT8(distance,dhDensity{dh,ii}.*1000,...
                sqrt(dhDensityVar{dh,ii}).*1000,1,{'Color',[0,0,0],'linewidth',1.5});
%                 sqrt(dhDensityVar{dh,ii}).*1000,1,{'Color',[0.5,0,0],'linewidth',1.5});
            hold on;
        end
        for rh = 1:nReflectionHorizon
            shadedErrorBarT8(distance,Density{rh,ii}.*1000,...
                sqrt(DensityVar{rh,ii}).*1000,1,{'Color',[0,0,0],'linewidth',1.5});
%                 sqrt(DensityVar{rh,ii}).*1000,1,{'Color',[0.5,0,0],'linewidth',1.5});
            hold on;
        end
        freezeColors
        axis ij
        axis tight
        grid on
        set(gca,'xticklabel',[])
        set(gca,'ytick',[350,375,400])
        title('Average Density (kg/m^{3})')
        subplot(3,1,3)
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
        set(gca,'ytick',[0.24,0.27,0.3,0.33])
        title('Average Annual Accumulation (m w.e. a^{-1})')
        set(gca,'xtick',1000.*[0, 20, 40, 60, 78])
        set(gca,'xticklabel',[0,20,40,60,78])
        xlabel('Distance (km)')
        set(findobj(gcf,'type','axes'),'FontName','FreeSerif','FontSize',12,...
            'FontWeight','Bold', 'LineWidth', 1);
    end
end
%% Compare Average Accumulation
if (isPickDepthHorizons || isLoadDepthHorizons) && (isPickAgeHorizons || isLoadIRH)
        for ii = 1:nFiles
        distance = Traverse{ii}./1000;
        figure(); hold on;
%          shadedErrorBarT8(distance,AverageAccumulation2{ii},...
%              sqrt(varAccumulation2{ii}),1,{'Color',[0,0,0],'linewidth',1.5});
%                  freezeColors

         shadedErrorBarT8(distance,AverageAccumulation3{ii},...
             sqrt(varAccumulation3{ii}),1,{'Color',[0,0,0],'linewidth',3});
                 freezeColors

                  for rh = 1:nReflectionHorizon
            shadedErrorBarT8(distance,AverageAccumulation{ii},...
                sqrt(SnowWaterEqvVar{rh,ii}),1,{'Color',[0.5,0,0],'linewidth',2});
                  end
        freezeColors
        axis tight
        axis ij
        grid on
        alpha(0.5)

%          title('Average Annual Accumulation (m w.e. a^{-1})')
         title({['\fontsize{14}Average Annual Accumulation'];...
             ['\fontsize{12}\color[rgb]{0.5 0 0}Jan.2015-Jan.2017 \color[rgb]{0 0 0}Jan.1984-Jan.2017']})
        set(gca,'xtick',[0, 20, 40, 60, 78])
        set(gca,'xticklabel',[0,20,40,60,78])
        set(gca,'fontsize',12,'fontweight','bold','layer','top')
        xlabel('Distance (km)')
        ylabel('(m w.e. a^{-1})')
        set(findobj(gcf,'type','axes'),'FontName','FreeSerif','FontSize',12,...
            'FontWeight','Bold', 'LineWidth', 1);
        
        % Compage GTC15 Density
          for ii = 1:nFiles
        distance = Traverse{ii}./1000;
        figure(); hold on;
%          shadedErrorBarT8(distance,AverageAccumulation2{ii},...
%              sqrt(varAccumulation2{ii}),1,{'Color',[0,0,0],'linewidth',1.5});
%                  freezeColors

         shadedErrorBarT8(distance,AverageAccumulation3{ii},...
             sqrt(SnowWaterEqvVar{1,ii}),1,{'Color',[0,0,0],'linewidth',2});

%              sqrt(varAccumulation3{ii}),1,{'Color',[0,0,0],'linewidth',3});

                 freezeColors

                  for rh = 1:nReflectionHorizon
            shadedErrorBarT8(distance,GTCaccumulation,...
                sqrt(SnowWaterEqvVar{rh,ii}),1,{'Color',[0.5,0,0],'linewidth',2});
            alpha(0.5)
                  end
        freezeColors
        axis tight
        axis ij
        grid on
%          title('Average Annual Accumulation (m w.e. a^{-1})')
         title({['\fontsize{14}Average Surface Mass Balance'];...
             ['\fontsize{12}\color[rgb]{0.5 0 0}GTC15 Density \color[rgb]{0 0 0}MXHL Density']})
        set(gca,'xtick',[0, 20, 40, 60, 78])
        set(gca,'xticklabel',[0,20,40,60,78])
        set(gca,'fontsize',12,'fontweight','bold')
        xlabel('Distance (km)')
        ylabel('(m w.e. a^{-1})')
          end
        
        figure();
        subplot(3,1,1)

        for rh = 1:nReflectionHorizon
            shadedErrorBarT8(distance,AverageAccumulation{ii},...
                sqrt(SnowWaterEqvVar{rh,ii}),1,{'Color',[0,0,0],'linewidth',1.5});
            hold on;
        end
        freezeColors
        axis tight
        axis ij
        grid on
        set(gca,'ytick',[0.24,0.27,0.3,0.33])
        title('Average Annual Accumulation (m w.e. a^{-1})')
        set(gca,'xtick',1000.*[0, 20, 40, 60, 78])
        set(gca,'xticklabel',[0,20,40,60,78])

        subplot(3,1,2)
            shadedErrorBarT8(distance,AverageAccumulation2{ii},...
                sqrt(varAccumulation2{ii}),1,{'Color',[0,0,0],'linewidth',1.5});            
            hold on;
        freezeColors
        axis tight
        axis ij
        grid on
        set(gca,'ytick',[0.24,0.27,0.3,0.33])
        title('Average Annual Accumulation (m w.e. a^{-1})')
        set(gca,'xtick',1000.*[0, 20, 40, 60, 78])
        set(gca,'xticklabel',[0,20,40,60,78])
        subplot(3,1,3)
            shadedErrorBarT8(distance,AverageAccumulation{ii},...
                sqrt(varAccumulation3{ii}),1,{'Color',[0,0,0],'linewidth',1.5});            
            hold on;
        freezeColors
        axis tight
        axis ij
        grid on
        set(gca,'ytick',[0.24,0.27,0.3,0.33])
        title('Average Annual Accumulation (m w.e. a^{-1})')
        set(gca,'xtick',1000.*[0, 20, 40, 60, 78])
        set(gca,'xticklabel',[0,20,40,60,78])
        xlabel('Distance (km)')
        set(findobj(gcf,'type','axes'),'FontName','FreeSerif','FontSize',12,...
            'FontWeight','Bold', 'LineWidth', 1);
    end
end