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
        ageStak = AgeModel{ii};
        DepositionAxis{ii} = (linspace(min(ageStak(:)),max(ageStak(:)),size(RadarDeposition{ii},1)))';
        DepositAxe = DepositionAxis{ii};
        % Grid Search for Depths Occurs in DepthtoDepositionConversion.m
        % Convert Deposition Time Image to Depth Image
        parfor (kk = 1:size(RadarDeposition{ii},2), nWorkers)
            % [Age,Depth] - Array to be resampled
            azcurve = [DepositAxe(:),zStak];
            [RadDepth(:,kk)] = timeDepthConversion(RadDeposit(:,kk),azcurve,depths(:,kk));

        end
        RadarDepth{ii} = RadDepth;
    end
    clear('RadarDeposition','RadDeposit','RadDepth','DepositionAxis','DepositAxe','zStak','ageStak','depths');
