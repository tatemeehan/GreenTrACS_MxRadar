%% Depth to Deposition
%
    % Allocation
    RadarDeposition = RadarDepth;  
    depositionAgeModel{ii} = cell(1,nFiles);
    DepositionAxis = cell(1,nFiles);
    for ii = 1:nFiles
        RadDeposition = zeros(size(RadarDeposition{ii}));
%         newAge = zeros(size(AgeModel{ii}));
        RadDepth = RadarDeposition{ii};
        % Axis for Deposition Image
        zStak = DepthAxis{ii};
        ageStak = AgeModel{ii};
        DepositionAxis{ii} = (linspace(min(ageStak(:)),max(ageStak(:)),size(RadarDeposition{ii},1)))';
        DepositAxe = DepositionAxis{ii};
        n = length(DepositAxe);
        depths = zeros(size(RadarDeposition{ii}));
        % Find the Max Depths
        % Grid Search for Depths of Isochrones at an evenly sampled rate.
        DepositAxeMat = repmat(DepositAxe,1,length(DepositAxe));
%         for (kk = 1:size(RadarDeposition{ii},2))
        parfor (kk = 1:size(RadarDeposition{ii},2),nWorkers)
            [~,tmpIx] = min(abs(DepositAxeMat - ageStak(:,kk)'));
            maxDepth = max(zStak(tmpIx));
            depths(:,kk) = linspace(0,maxDepth,n);
        end
        % Do the Converison
        parfor (kk = 1:size(RadarDeposition{ii},2), nWorkers)
            % [Depth, Age] - Array to be resampled
            zacurve = [depths(:,kk),ageStak(:,kk)];
            % Perform Depth Conversion
            [RadDeposition(:,kk)] = timeDepthConversion(RadDepth(:,kk),zacurve,DepositAxe);
            % Convert Age-Depth Model 
            % Yes this Works, but it is simpler to compute this directly
%             [newAge(:,kk)] = timeDepthConversion(ageStak(:,kk),zacurve,DepositAxe);
        end
        RadarDeposition{ii} = RadDeposition;
%         depositionAgeModel{ii} = newAge;
        depositionAgeModel{ii} = ones(1,size(RadarDeposition{ii},2)).*DepositionAxis{ii};
    end
    clear('RadDeposition','RadDepth','DepositAxe','zStak','ageStak');
