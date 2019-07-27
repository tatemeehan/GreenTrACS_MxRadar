%% Caclulate Stratigraphic Age Residual
% The residual of the picked age-horizons is measured optionally against 
% the nearby GreenTrACS firn core or the starting Age-Depth model.

% Determine Number of Age Horizons
nAgeHorizon = size(isochronePick,2);

% Calculate Residual
ageResidual = isochronePick;
datumAge = zeros(nAgeHorizon,nFiles);
warning('off','MATLAB:mir_warning_changing_try_catch');
try iceCoreIx;
catch
    iceCoreIx = 1:100;
end

% Cacluate the Age Residual each row of the matrix should have the same age
for ii = 1:nFiles
    dAge = mean(diff(DepositionAxis{ii}));
    for jj = 1;% chan (Loop over chan for Mx processing)
             % Test for Age samples or Index
            if all(max([isochronePick{jj,:,ii}]) > 25) % 25 samples/age threshold
                isConvert = 1;
            else
                isConvert = 0;
            end
        for ah = 1:nAgeHorizon
            if isConvert
            % Convert Samples to Age
            isochronePick{jj,ah,ii} = isochronePick{jj,ah,ii}.*dAge;
            end
            % Determine Measurement Datum
            if isGreenTracsFirnCore
                % ageDepth : Age is the linear Axis!
                GTCageDepth = interp1(abs(GTCaz(:,1)-GTCaz(1,1)),GTCaz(:,2),DepositionAxis{ii},'pchip');
                GTCageDepth = [DepositionAxis{ii},GTCageDepth];
                % Grab the mean Starting Age-Depth Model at iceCoreIx
                % Interpolate these depths to the Deposition Axis
                % With this new axis find the Depths of the picks
                %Find equivalent depths on GTC ageDepth axis
                % Compute residual ages of equivalent depths
                tmpA = mean(AgeModel{ii}(:,iceCoreIx),2);
                tmpZ = mean(DepthMatrix{ii}(:,iceCoreIx),2);
                tmp = interp1(tmpA,tmpZ,DepositionAxis{ii},'pchip');
                modAgeDepth = [DepositionAxis{ii},tmp];
                isochroneDepth(ah,ii) = modAgeDepth(round(mean(isochronePick{jj,ah,ii}(iceCoreIx)./dAge)),2);
                [~,datumIx] = min(abs(GTCageDepth(:,2) - isochroneDepth(ah,ii)));
            else
                [~,datumIx] = min(abs(mean(isochronePick{jj,ah,ii}(iceCoreIx))-DepositionAxis{ii}));
            end
            datumAge(ah,ii) = DepositionAxis{ii}(datumIx);
            % Calculate Misfit
            ageResidual{jj,ah,ii} = isochronePick{jj,ah,ii} - datumAge(ah,ii);
        end
    end
    clear ('tmp','tmpA','tmZ','modAgeDepth','datumIx','isConvert');%,'isochroneDepth','GTCageDepth',
end