%% Visualize Data Gathers
% This Routine allows the user to display the GPR Sections. 

        if isPlotOffsetGathers
            % AGC Gain for Plotter
            for ii = 1:nFiles
                for jj = chan
                    plotRadar{jj,ii} = AGCgain(Radar{jj,ii},50,2);
                end
            end

            if isTopMute
            % Top Mute for Plotter
            for ii = 1:nFiles
                for jj = chan
                    plotRadar{jj,ii}(1:xcorrWindow(2,jj),:) = ones(xcorrWindow(2,jj),1)*Radar{jj,ii}(1,:);
                end
            end
            end
            
            % Plot Travel-Time Picks
%             if isLoadTimeHorizons
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
