%% Interpolate Residual Update Model
% The residual of the picked isochrone horizons is interpolated by 1D shape 
% preserving polynomials in within a loop, gaussian smoothing is applied.
% The 2D residual model is summed to the prior age model as a perturbation
% update. 
%
deltaAgeModel = cell(1,nFiles);
depositionAgeModel = cell(1,nFiles);

    for ii = 1:nFiles
        depositionAgeModel{ii} = ones(1,size(RadarDeposition{ii},2)).*DepositionAxis{ii};
        isochroneResidual = [ageResidual{1,:,ii}]';
        npicks = size(isochroneResidual,2);
        % Append Surface Data ~ Assume no misfit at present day
        isochroneResidual = [zeros(1,npicks);isochroneResidual];
        tmpDatum = [0;datumAge(:,ii)];
        % 1D Interpolation
        for jj = 1:npicks
         % Depth-Age Model Update
         deltaAgeModel{ii}(:,jj) = interp1(tmpDatum,isochroneResidual(:,jj),DepositionAxis{ii},'pchip',0);
        end
        % Smooth deltaAgeModel
        deltaAgeModel{ii} = imgaussfilt(deltaAgeModel{ii},[50,5]);
        depositionAgeModel{ii} = depositionAgeModel{ii} - deltaAgeModel{ii};
    end
    clear('isochroneResidual','npicks','tmpDatum')
    
    
    %% Convert Radar Traces to Updated Age/Isochrone Model
    % The age/deposition radargram is converted by the updated age-depth 
    % model to perform residual trace flattening. The age-depth model is
    % updated by the residual corrections and the average annual
    % accumulation is re-calculated.
    
    for ii = 1:nFiles
        RadDeposition = zeros(size(RadarDeposition{ii}));
        RadDeposit = RadarDeposition{ii};
        tStak = tStack{ii};
        zStak = zStack{ii};
        % Axis for Deposition Image
        ageStak = pseudoAgeModel{ii};
        DepositAxe = DepositionAxis{ii};
        DepthAxe = DepthAxis{ii};
        n = length(DepositAxe);
        ageUpdate = depositionAgeModel{ii};
        ages = ones(1,size(RadarDeposition{ii},2)).*DepositionAxis{ii};%zeros(size(RadarDeposition{ii}));
        newAge = zeros(size(zStak));
        % Do the Converison
        parfor (kk = 1:size(RadarDeposition{ii},2), nWorkers)
%         for kk = 1:size(RadarDeposition{ii},2)
            % Perform Trace Flattening
            RadDeposition(:,kk) = interp1(ageUpdate(:,kk),RadDeposit(:,kk),ages(:,kk),'linear');
            % Update Age-Depth Model
            % Create the Conversion Axis Time 2 Age (TA)
            tmpAxTA = interp1(ageStak(:,kk),tStak(:,kk),DepositAxe,'linear');
            nanIx = isnan(tmpAxTA);
            tmpAxTA = tmpAxTA(~nanIx);
            % Create the Conversion Axis Age 2 Depth (TZ)
            tmpAxTZ = interp1(tStak(:,kk),zStak(:,kk),tmpAxTA,'linear');
            mData = length(tmpAxTZ);
            % Use interp1 Alg Here
            newAge(:,kk) = interp1(tmpAxTZ,ageUpdate(1:mData,kk),DepthAxe);
        end
        % Fill Interpolation NaNs
        nanIx = isnan(RadDeposition);
        RadDeposition(nanIx) = 0;
        RadarDeposition{ii} = RadDeposition;
        newAge = abs(imgaussfilt(inpaint_nans(newAge,1),[1,50]));
        newAge(1,:) = 0;
        updateAgeModel{ii} = newAge;
        dUpdateAgeModel = [eps.*ones(1,size(newAge,2));diff(newAge)];
%                 dUpdateAgeModel = [eps.*ones(1,size(updateAgeModel{ii},2));diff(updateAgeModel{ii})];

        updatePseudoAgeModel{ii} = (StackingVelocityModel{ii}.*(updateAgeModel{ii})./DepthMatrix{ii}).*(tStack{ii}./2);
        updatePseudoAgeModel{ii}(1,:) = 0; % Set NaNs to Zero
%         tmpDiff = diff(DepthMatrix{ii});
%         dUpdateAgeModel = [tmpDiff(1,:);tmpDiff];
        
        %%% Causal Integration Inversion for Instantaneous Accumulation %%% 
        %AverageAccumulation{ii} = mean((DepthMatrix{ii}(accumIx:end,:).*AvgDensityModel{ii}(accumIx:end,:))./(updateAgeModel{ii}(accumIx:end,:)+eps));        
        AverageAccumulationMatrix = (DepthMatrix{ii}.*AvgDensityModel{ii})./(updateAgeModel{ii}+eps);
        instantAccum = NaN(size(updateAgeModel{ii}));
        A1 = 2017; % Upper Year Bound
        % Lower Year Bound A2 = 1991
        A2 = round(abs(max(datumAge)-(str2num(Year{1})+dayofyear/365)));
        for kk = 1:size(updateAgeModel{ii},2)
            % Find indicies for Age Range A1 - A2
            ix = find(abs(updateAgeModel{ii}(:,kk)-(str2num(Year{1})+dayofyear/365)) < A1 ...
                & abs(updateAgeModel{ii}(:,kk)-(str2num(Year{1})+dayofyear/365)) >= A2);
            % Compute Weights for Accumulation
            W = dUpdateAgeModel(ix,kk)'.*tril(ones(length(ix)));
            Gw = W./sum(W,2);
            dw = AverageAccumulationMatrix(ix,kk);
            % Compute the Derivative
            m = Gw\dw;
            instantAccum(ix,kk) = m;
        end
        AverageAccumulation2{ii} = nanmean(instantAccum);
        varAccumulation2{ii} = nanvar(instantAccum);
        initAccum = repmat(AverageAccumulation{ii},size(updateAgeModel{ii},1),1);
        instantDoom = instantAccum;
        for kk = 1 :length(AverageAccumulation{ii})
            nanIx = isnan(instantAccum(:,kk));
            instantDoom(nanIx,kk) = AverageAccumulation{ii}(kk);
        end
    end
    clear('n','RadDeposition','RadDeposit','DepositAxe','DepositAxeMat',...
        'ageStak','aacurve','ageUpdate','azcurve','newAge','accumIx',...
        'dUpdateAgeModel','W','Gw','dw','m','AverageAccumulationMatrix',...
        'instantAccum');


    
    