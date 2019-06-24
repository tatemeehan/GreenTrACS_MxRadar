    %% Convert Radar Traces to Updated Age/Isochrone Model
    % The age/deposition radargram is converted by the updated age-depth 
    % model to perform residual trace flattening.
    
    for ii = 1:nFiles
        RadDeposition = zeros(size(RadarDeposition{ii}));
        RadDeposit = RadarDeposition{ii};
        % Axis for Deposition Image
        ageStak = AgeModel{ii};
        DepositAxe = DepositionAxis{ii};
        ageUpdate = depositionAgeModel{ii};
        % Do the Converison
        parfor (kk = 1:size(RadarDeposition{ii},2), nWorkers)
%         for (kk = 1:size(updateRadarDeposition{ii},2))
            % [Age, updateAge] - Array to be resampled
            aacurve = [DepositAxe,ageUpdate(:,ii)];
            % Perform Trace Unflattening
            [RadDeposition(:,kk)] = timeDepthConversion(RadDeposit(:,kk),aacurve,ageUpdate(:,kk));
        end
        RadarDeposition{ii} = RadDeposition; 
    end