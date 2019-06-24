%% Interpolate Residual Update Model
% The residual of the picked isochrone horizons is interpolated by 1D shape 
% preserving polynomials in within a loop, gaussian smoothing is applied.
% The 2D residual model is summed to the prior age model as a perturbation
% update. 
%
deltaAgeModel = cell(1,nFiles);
    for ii = 1:nFiles
        isochroneResidual = [ageResidual{1,:,ii}]';
        npicks = size(isochroneResidual,2);
        % Append Surface Data ~ Assume no misfit at present day
        isochroneResidual = [zeros(1,npicks);isochroneResidual];
        tmpDatum = [0;datumAge(:,ii)];
        % 1D Interpolation
        for jj = 1:npicks
         % Depth-Age Model Update
         deltaAgeModel{ii}(:,jj) = interp1(tmpDatum,isochroneResidual(:,jj),DepositionAxis{ii},'pchip');
        end
        % Smooth deltaAgeModel
        deltaAgeModel{ii} = imgaussfilt(deltaAgeModel{ii},50);
        depositionAgeModel{ii} = depositionAgeModel{ii} - deltaAgeModel{ii};
    end
    clear('isochroneResidual','npicks','tmpDatum')
    
    
    %% Convert Radar Traces to Updated Age/Isochrone Model
    % The age/deposition radargram is converted by the updated age-depth 
    % model to perform residual trace flattening.
    
    for ii = 1:nFiles
        RadDeposition = zeros(size(RadarDeposition{ii}));
        RadDeposit = RadarDeposition{ii};
        % Axis for Deposition Image
        ageStak = AgeModel{ii};
        DepositAxe = DepositionAxis{ii};
        DepthAxe = DepthAxis{ii};
        n = length(DepositAxe);
        ageUpdate = depositionAgeModel{ii};
        ages = zeros(size(RadarDeposition{ii}));
        newAge = ages;
        % Find the Max Ages use linear sampling
        % Grid Search for Ages of Isochrones at an evenly sampled rate.
        DepositAxeMat = repmat(DepositAxe,1,length(DepositAxe));
%         for (kk = 1:size(RadarDeposition{ii},2))
        parfor (kk = 1:size(RadarDeposition{ii},2),nWorkers)
            [~,tmpIx] = min(abs(DepositAxeMat - ageStak(:,kk)'));
            maxAge = max(DepositAxe(tmpIx));
            ages(:,kk) = linspace(0,maxAge,n);
        end
        % Do the Converison
        parfor (kk = 1:size(RadarDeposition{ii},2), nWorkers)
%         for (kk = 1:size(updateRadarDeposition{ii},2))
            % [Age, updateAge] - Array to be resampled
            aacurve = [ages(:,kk),ageUpdate(:,kk),];
            % Perform Trace Flattening
            [RadDeposition(:,kk)] = timeDepthConversion(RadDeposit(:,kk),aacurve,DepositAxe);
            % Update Age-Depth Model            
            azcurve = [DepositAxe,DepthAxe];
            [newAge(:,kk)] = timeDepthConversion(ageUpdate(:,kk),azcurve,depths(:,kk));
        end
        RadarDeposition{ii} = RadDeposition; 
        updateAgeModel{ii} = newAge;
    end
    clear('n','RadDeposition','RadDeposit','DepositAxe','DepositAxeMat','ageStak','aacurve','ageUpdate','azcurve','newAge');


    
    