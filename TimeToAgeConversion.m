%% Time To Stratigraphic Age Conversion
    % Allocation
    RadarDeposition = RadarStack;  
    depositionAgeModel{ii} = cell(1,nFiles);
    DepositionAxis = cell(1,nFiles);
    pseudoAgeModel{ii} = cell(1,nFiles);

    for ii = 1:nFiles
        DepositionAxis{ii} = [linspace(min(AgeModel{ii}(:)),max(AgeModel{ii}(:)),length(AgeModel{ii}(:,1)))]';
        % Create Age-Time Model
        pseudoAgeModel{ii} = (StackingVelocityModel{ii}.*(AgeModel{ii})./DepthMatrix{ii}).*(tStack{ii}./2);
        pseudoAgeModel{ii}(1,:) = 0; % Set NaNs to Zero
        tmpModel = pseudoAgeModel{ii};

        [m,n] = size(RadarStack{ii});
        q = 2; % Quantity times the current number of samples
        RadDeposit = zeros(q.*m,n);
        RadStack = RadarStack{ii};
        DepositAxe = DepositionAxis{ii};
        % ReSample in Age Domain
        DepositAxe = [linspace(min(DepositAxe),max(DepositAxe),q.*length(DepositAxe))]';
        DepositionAxis{ii} = DepositAxe;
        zStak = zStack{ii};
        tStak = tStack{ii};
        parfor (kk = 1:size(RadarStack{ii},2), nWorkers)
%         for kk = 1:size(RadarStack{ii},2)
%             % [Time, Age] - Array to be resampled
%             tacurve = [tStak(:,kk),tmpModel(:,kk)];
%             % Perform Depth Conversion
%             [RadDeposit(:,kk)] = timeDepthConversion(RadStack(:,kk),tacurve,DepositAxe);
            % Create the Conversion Axis Time 2 Age (TA)
            tmpAxTA = interp1(tmpModel(:,kk),tStak(:,kk),DepositAxe,'linear');
            nanIx = isnan(tmpAxTA);
            tmpAxTA=tmpAxTA(~nanIx);
            % Use interp1 Alg Here
            tmpTrc = interp1(tStak(:,kk),RadStack(:,kk),tmpAxTA,'linear');
            % Pad Interpolation with Zeros
            RadDeposit(:,kk) = [tmpTrc;zeros(q.*m - length(tmpTrc),1)];
        end
        RadarDeposition{ii} = RadDeposit;
    end
    clear('RadDepth','RadStack','DepthAxe','zStak','tStak','tmpModel','RadDeposit','DepositAxe');