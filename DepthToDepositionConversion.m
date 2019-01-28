%% Depth to Deposition
%
    % Allocation
    RadarDeposition = RadarDepth;  
    DepositionAxis = cell(1,nFiles);
    for ii = 1:nFiles
        RadDeposition = zeros(size(RadarDeposition{ii}));
        RadDepth = RadarDeposition{ii};
        % Axis for Deposition Image
        zStak = DepthAxis{ii};
        ageStak = AgeModel{ii};
        DepositionAxis{ii} = (linspace(min(ageStak(:)),max(ageStak(:)),size(RadarDeposition{ii},1)))';
        DepositAxe = DepositionAxis{ii};
        n = length(DepositAxe);
        depths = zeros(size(RadarDeposition{ii}));
        % Find the Max Depths
        % Grid Search for Depths or Isochrones at an evenly sampled rate.
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
%             zacurve = [zStak(:),ageStak(:,kk)];
            % Perform Depth Conversion
            [RadDeposition(:,kk)] = timeDepthConversion(RadDepth(:,kk),zacurve,DepositAxe);
        end
        RadarDeposition{ii} = RadDeposition;
    end
    clear('RadDeposition','RadDepth','DepositAxe','zStak','ageStak');
