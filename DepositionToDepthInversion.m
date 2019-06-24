%% Deposition to Depth Inversion
%
% Allocation
RadarDepth = RadarDeposition;
% Change to Age Depth Matrix
    for ii = 1:nFiles
        RadDeposit = RadarDeposition{ii};
        RadDepth = zeros(size(RadarDeposition{ii}));
        % Axis for Deposition Image
        zStak = DepthAxis{ii};
        if isPickAgeHorizons || isLoadIRH
            % The Deposition Age Model needs conversion to Age-Depth!!!
            ageStak = depositionAgeModel{ii};
            DepositionAxis{ii} = (linspace(0,max(ageStak(:)),size(RadarDeposition{ii},1)))';
            DepositAxe = DepositionAxis{ii};
%             depths = zeros(size(RadarDeposition{ii}));
%             % Find the Max Depths
%             % Grid Search for Depths of Isochrones at an evenly sampled rate.
%             DepositAxeMat = repmat(DepositAxe,1,length(DepositAxe));
%             %         for (kk = 1:size(RadarDeposition{ii},2))
%             parfor (kk = 1:size(RadarDeposition{ii},2),nWorkers)
%                 [~,tmpIx] = min(abs(DepositAxeMat - ageStak(:,kk)'));
%                 maxDepth = max(zStak(tmpIx));
%                 depths(:,kk) = linspace(0,maxDepth,n);
%             end
        else
        ageStak = AgeModel{ii};
        DepositionAxis{ii} = (linspace(0,max(ageStak(:)),size(RadarDeposition{ii},1)))';
        DepositAxe = DepositionAxis{ii};
        end
        % Grid Search for Depths Occurs in DepthtoDepositionConversion.m
        % Convert Deposition Time Image to Depth Image
        parfor (kk = 1:size(RadarDeposition{ii},2), nWorkers)
            % [Age,Depth] - Array to be resampled
            azcurve = [DepositAxe(:),zStak];
            [RadDepth(:,kk)] = timeDepthConversion(RadDeposit(:,kk),azcurve,depths(:,kk));

        end
        RadarDepth{ii} = RadDepth;
    end
    clear('RadDeposit','RadDepth','DepositAxe','zStak','ageStak');%'RadarDeposition','DepositionAxis','depths'
