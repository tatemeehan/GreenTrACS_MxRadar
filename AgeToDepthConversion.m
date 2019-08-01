%% Stratigraphic Age to Depth Conversion
%
% Allocation
RadarDepth = RadarDeposition;
% Change to Age Depth Matrix
    for ii = 1:nFiles
        RadDeposit = RadarDeposition{ii};
        [m,n] = size(RadarDeposition{ii});
        m = m/q;
        RadDepth = zeros(m,n);
        DepositAxe = DepositionAxis{ii};
        % ReSample in Age Domain
%         DepositAxe = [linspace(min(DepositAxe),max(DepositAxe),q.*length(DepositAxe))]';
        zStak = zStack{ii};
        tStak = tStack{ii};
%         tmpModel = pseudoAgeModel{ii};
        tmpModel = updatePseudoAgeModel{ii};
        % Interpolation Axes
        AxZ = linspace(min(zStak(:)),max(zStak(:)),m);
        % Axis for Deposition Image
%         if isPickAgeHorizons || isLoadIRH
%             ageStak = depositionAgeModel{ii};
%             DepositionAxis{ii} = (linspace(0,max(ageStak(:)),size(RadarDeposition{ii},1)))';
%             DepositAxe = DepositionAxis{ii};
%         else
%         end
        % Convert Stratigraphic Age Image to Depth Image
        parfor (kk = 1:size(RadarDeposition{ii},2), nWorkers)
%         for kk = 1:size(RadarDeposition{ii},2)
            % Create the Conversion Axis Time 2 Age (TA)
            tmpAxTA = interp1(tmpModel(:,kk),tStak(:,kk),DepositAxe,'linear');
            nanIx = isnan(tmpAxTA);
            tmpAxTA = tmpAxTA(~nanIx);
            % Create the Conversion Axis Age 2 Depth (TZ)
            tmpAxTZ = interp1(tStak(:,kk),zStak(:,kk),tmpAxTA,'linear');
            mData = length(tmpAxTZ);

            % Use interp1 Alg Here
            tmpTrc = interp1(tmpAxTZ,RadDeposit(1:mData,kk),AxZ);
            % Pad Interpolation with Zeros
            RadDepth(:,kk) = [tmpTrc;zeros(m - length(tmpTrc),1)];

            % [Age,Depth] - Array to be resampled
%             azcurve = [tmpModel(:,kk),zStak(:,kk)];
%             [RadDepth(:,kk)] = timeDepthConversion(RadDeposit(:,kk),azcurve,DepthAxe);
        end
        nanIx = isnan(RadDepth);
        RadDepth(nanIx) = 0;
        RadarDepth{ii} = RadDepth;
    end
    clear('RadDeposit','RadDepth','DepositAxe','zStak','ageStak','tmpModel','tStak','AxZ');%'RadarDeposition','DepositionAxis','depths'