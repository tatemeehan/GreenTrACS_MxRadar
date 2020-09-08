%% Caclulate Stratigraphic Depth Residual
% The residual of the picked depth-horizons is measured optionally against 
% the nearby GreenTrACS firn core.
% depthPick(21) = [];
% depthPick(1) = [];
% Determine Number of Age Horizons
nZHorizon = size(depthPick,2);

% Calculate Residual
isoResidual = depthPick;
agePick = depthPick;%zeros(nZHorizon,nFiles);
warning('off','MATLAB:mir_warning_changing_try_catch');
try iceCoreIx;
catch
    iceCoreIx = 1:100;
end

% Cacluate the Age Residual each row of the matrix should have the same age
for ii = 1:nFiles
    dZ = mean(diff(DepthAxis{ii}));
    for jj = 1;% chan (Loop over chan for Mx processing)
             % Test for Age samples or Index
            if all(max([depthPick{jj,:,ii}]) > 20) % 20 samples/depth threshold
                isConvert = 1;
            else
                isConvert = 0;
            end
        for zh = 1:nZHorizon
            if isConvert
            % Convert Samples to Age
            depthPick{jj,zh,ii} = depthPick{jj,zh,ii}.*dZ;
            end
            % Determine Measurement Datum
            if isGreenTracsFirnCore
                % depthAge : Depth is the linear Axis!
                GTCdepthAge = interp1(abs(GTCdepthAge(:,1)-GTCdepthAge(1,1)),GTCdepthAge(:,2),DepthAxis{ii},'pchip');
                GTCdepthAge = [DepthAxis{ii},GTCdepthAge];
                % Average Core Location Isochrone Depth as Datum
                isochroneDepth(zh,ii) = mean(depthPick{jj,zh,ii}(iceCoreIx));
                [~,datumIx] = min(abs(GTCdepthAge(:,1) - isochroneDepth(zh,ii)));
            else
                [~,datumIx] = min(abs(mean(depthPick{jj,zh,ii}(iceCoreIx))-DepositionAxis{ii}));
            end
            % Grab the mean updated Age-Depth Model at iceCoreIx
%             ageDatum(zh,ii) = mean(updateAgeModel{ii}(datumIx,iceCoreIx));
               tmpIx = sub2ind(size(updateAgeModel{ii}),round(depthPick{jj,zh,ii}./dZ),...
                   [1:length(depthPick{jj,zh,ii})]');
            agePick{jj,zh,ii} = updateAgeModel{ii}(tmpIx);
            % Calculate Misfit -- Subtract Residuals from Model for Update! 
            datumDepthAge(zh) =  abs(GTCdepthAge(datumIx,2)-GTCdepthAge(1,2));
            isoResidual{jj,zh,ii} = agePick{jj,zh,ii} - abs(GTCdepthAge(datumIx,2)-GTCdepthAge(1,2));
        end
        % Interpolate Perterbations to Depth
        tmpResidual = [isoResidual{:}];tmpResidual = [zeros(length(tmpResidual),1),tmpResidual];
        tmpDepthPicks = [depthPick{:}];tmpDepthPicks = [zeros(length(tmpResidual),1),tmpDepthPicks];
        tmpAxis = DepthAxis{ii};
        perturbations = zeros(size(updateAgeModel{ii}));
        tmpAgeModel = perturbations;
        parfor (kk = 1:length(updateAgeModel{ii}), nWorkers)
%         for kk = 1:length(updateAgeModel{ii})
        perturbations(:,kk) = interp1(tmpDepthPicks(kk,:),tmpResidual(kk,:),tmpAxis,'linear','extrap');
        end
        % Smooth Perturbations
        perturbations = imgaussfilt(perturbations,[0.1,50]);
        for kk = 1:length(updateAgeModel{ii})
        % Re-Update AgeDepth Model
        tmpAgeModel(:,kk) = updateAgeModel{ii}(:,kk) - perturbations(:,kk);
        % Correction to Surface Gradients
        tmpIx = find(diff(sign(tmpAgeModel(:,kk))),1,'last');
        if tmpIx > 1
        tmpIx = (20.*tmpIx) + 1;
        surfIx = 1:tmpIx;
        tmpAgeModel(surfIx,kk) = linspace(0,tmpAgeModel(tmpIx,kk),tmpIx);
        end
        end
        bestAgeModel{ii} = tmpAgeModel;
        % Update Age Picks
        for zh = 1:nZHorizon
        % Grab the mean updated Age-Depth Model at iceCoreIx
%             ageDatum(zh,ii) = mean(updateAgeModel{ii}(datumIx,iceCoreIx));
               tmpIx = sub2ind(size(bestAgeModel{ii}),round(depthPick{jj,zh,ii}./dZ),...
                   [1:length(depthPick{jj,zh,ii})]');
            agePick{jj,zh,ii} = bestAgeModel{ii}(tmpIx);
        end
        clear ('tmpResidual','tmpDepthPicks','tmpAxis','tmpIx','tmpAgeModel')
        % 1st order Finite Difference 
        tmpDiffA = diff(bestAgeModel{ii});  % da
        tmpDiff = diff(DepthMatrix{ii});    % dz
        dUpdateAgeModel = [tmpDiff(1,:);tmpDiff]; % Colocation
%         % Compute the Integrated Average Accumulation
%         AverageAccumulationMatrix = (DepthMatrix{ii}.*AvgDensityModel{ii})./(bestAgeModel{ii}+eps);
%         % Compute dz/da for instantaneous accumulation rate
        GTCdensity = interp1(GTCzp(:,1),GTCzp(:,2),DepthAxis{ii},'pchip');
        GTCdensityModel = repmat(GTCdensity,1,n);
        GTCaccum = (dUpdateAgeModel.*GTCdensityModel)./([tmpDiffA(1,:);tmpDiffA]);
        instantAccumO1 = (dUpdateAgeModel.*DensityModel{ii})./([tmpDiffA(1,:);tmpDiffA]);
% 2nd Order Finite Difference
tic
instantSMB = zeros(length(DepthMatrix{ii}(:,1)),length(updateAgeModel{ii}));
parfor (kk = 1:length(updateAgeModel{ii}), nWorkers)
% for kk = 1:length(updateAgeModel{ii})
    tmp = zeros(length(DepthMatrix{ii}(:,1)),1);
    
for ll = 1:length(DepthMatrix{ii}(:,1))
    locZ = DepthMatrix{ii}(ll,kk);
    [~,ptsIx] = mink(abs(locZ-DepthMatrix{ii}(:,kk)),3);%find(abs(locZ-DepthMatrix{ii}(:,kk))>=0,3);
    ptsIx = sort(ptsIx);
    ptsZ = DepthMatrix{ii}(ptsIx,kk);
%     cz = fdweights(locZ,ptsZ,1);
%     dz(ll) = cz*ptsZ;
    locA = bestAgeModel{ii}(ll,kk);
    ptsA = bestAgeModel{ii}(ptsIx,kk);
    ca = fdweights(locA,ptsA,1);
%     da(ll) = ca*ptsA;
    tmp(ll) = (ca(2,:)*ptsZ)*DensityModel{ii}(ll,kk);
end
instantSMB(:,kk) = tmp;
%% 
end

        % Despike Mathy Noise
        cutoff = quantile(instantSMB(:),[.005,.995]);
        tmpIx = find(instantSMB<cutoff(1) | instantSMB > cutoff(2));
        instantSMB(tmpIx) = NaN;
%         testInstant = inpaint_nans(testInstant,0); % Quite Slow
        instantSMB = fillmissing(instantSMB,'linear');
        clear('cutoff','tmpIx')
        A1 = 2017; % Upper Year Bound
        % Lower Year Bound A2 = 1984
        A2 = round(abs(max(datumDepthAge)-(str2num(Year{1})+dayofyear/365)));
        Annuals = A2:A1;
        tmpAccum = instantSMB;
        m = size(bestAgeModel{ii},1);
        n = size(bestAgeModel{ii},2);
        reportAccum = zeros(n,1);
        varAccum = zeros(n,1);
        bin = zeros(length(Annuals)-1,n);
        gtc = bin;binVar = bin; gtcVar = bin;binN = bin; gtcN= bin;
        for kk = 1:n
            % Find indicies for Age Range A1 - A2
            ix = find(abs(bestAgeModel{ii}(:,kk)-(str2num(Year{1})+dayofyear/365)) < A1 ...
                & abs(bestAgeModel{ii}(:,kk)-(str2num(Year{1})+dayofyear/365)) >= A2);
            tmpAgeModel = bestAgeModel{ii}(ix,kk);
            tmpA = instantSMB(ix,kk);
            tmpGTCa = GTCaccum(ix,kk);

            % Pad update to Instantaeous Accumulation Model with intial 
            tmpAccum(~ismember(1:m,ix),kk) = AverageAccumulation{ii}(kk);% tmpAccum is instantDoom
            for ll = 1:length(Annuals)-1
                tmpIx = find(tmpAgeModel<abs(Annuals(ll)-(str2num(Year{1})+dayofyear/365))&...
                    tmpAgeModel>=abs(Annuals(ll+1)-(str2num(Year{1})+dayofyear/365)));
                bin(ll,kk) = mean(tmpA(tmpIx));
                binVar(ll,kk) = var(tmpA(tmpIx));
                binN(ll,kk) = length(tmpA(tmpIx));
                gtc(ll,kk) = mean(tmpGTCa(tmpIx));
                gtcVar(ll,kk) = var(tmpGTCa(tmpIx));
                gtcN(ll,kk) = length(tmpGTCa(tmpIx));
            end
            
            % Monte Carlo Sampling for Decaldal Mean accumulation and var
            mcAccum = zeros(100,1);
            for ll = 1:100
                mcAccum(ll) =mean(datasample(bin(:,kk),10,'replace',false));
                mcGTCaccum(ll) = mean(datasample(gtc(:,kk),10,'replace',false));
            end
%             reportAccum(kk) = mean(instantAccum(ix,kk));
%             varAccum(kk) = var(instantAccum(ix,kk));
            reportAccum(kk) = mean(mcAccum);
            varAccum(kk) = var(mcAccum);
%             GTCaccumulation(kk) = mean(GTCaccum(ix,kk));
%             GTCvarAccumulation(kk) = var(GTCaccum(ix,kk));
            GTCaccumulation(kk) = mean(mcGTCaccum);
            GTCvarAccumulation(kk) = var(mcGTCaccum);
        end
        reportAccum = nonParametricSmooth(Traverse{ii},reportAccum,Traverse{ii},251);
        varAccum = nonParametricSmooth(Traverse{ii},varAccum,Traverse{ii},251);
        
        GTCaccumulation = nonParametricSmooth(Traverse{ii},GTCaccumulation,Traverse{ii},251);
        GTCvarAccumulation = nonParametricSmooth(Traverse{ii},GTCvarAccumulation,Traverse{ii},251);

        AverageAccumulation3{ii} = reportAccum;
        varAccumulation3{ii} = varAccum;
        % Calculate Experimental Confidence Region
        confidance = zeros(n,1);
        for kk = 1:n
            if AverageAccumulation3{ii}(kk) > AverageAccumulation{ii}(kk)+sqrt(SnowWaterEqvVar{1,ii}(kk))
                confidance(kk) = 1;
            elseif AverageAccumulation3{ii}(kk) < AverageAccumulation{ii}(kk)-sqrt(SnowWaterEqvVar{1,ii}(kk))
                confidance(kk) = 1;
            end
        end
            confidanceIx = find(confidance);
       CI = 1 - (length(confidanceIx)./length(AverageAccumulation3{ii}));
       clear('reportAccum','varAccum','m','n');
    end
    clear ('tmp','tmpA','tmZ','datumIx','isConvert','tmpIx','tmpResidual');%,'isochroneDepth','GTCageDepth',
end