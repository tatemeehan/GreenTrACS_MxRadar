%% Data Static Shift
% This Routine Applies the Necessary Trace Rotation to the Data Gathers
% that correct for the Digital Time-Sampling Errors on the rocording.

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
for ii = 1:nFiles
    parfor (jj = 1:nChan, nWorkers)
        % Remove AirWave 0verhead Time
        Overhead = Radar{jj,ii};
        OverheadIx = round(AirTo{ii,jj}./dt);
        Radar{jj,ii} = [Overhead(OverheadIx+1:end,:);Overhead(1:OverheadIx,:)];
        % Create Temporary Dimension Handle
        tmpLength = size(Radar{jj,ii},1);
        noise = mean(Radar{jj,ii}(1,:));
        directNoise = mean(directRadar{jj,ii}(1,:));
        for kk = 1:size(deltaT,2)
            
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