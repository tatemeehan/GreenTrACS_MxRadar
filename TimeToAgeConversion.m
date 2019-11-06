%% Time To Stratigraphic Age Conversion
    % Allocation
    RadarDeposition = RadarStack;  
    depositionAgeModel{ii} = cell(1,nFiles);
    DepositionAxis = cell(1,nFiles);
    pseudoAgeModel{ii} = cell(1,nFiles);

    for ii = 1:nFiles
        % Create Age-Time Model
%         pseudoAgeModel{ii} = (StackingVelocityModel{ii}.*(AgeModel{ii})./DepthMatrix{ii}).*(tStack{ii}./2);
%         pseudoAgeModel{ii} = (StackingVelocityModel{ii}.*(AgeModel{ii})./zStack{ii}).*(tStack{ii}./2);
%         pseudoAgeModel{ii}(1,:) = 0; % Set NaNs to Zero
%         tmpModel = pseudoAgeModel{ii};
        % Using Depth to Time Conversion
        ageMod = AgeModel{ii};
        zStak = zStack{ii};
        tStak = tStack{ii};
        zAxe = DepthMatrix{ii};
        ageT = zeros(size(ageMod));
        parfor (kk = 1:size(RadarStack{ii},2), nWorkers)
            % Conversion Axis Depth 2 Time (ZT) is zStacK
            % Use interp1 Alg Here
            tmpAgeT = interp1(zAxe(:,kk),ageMod(:,kk),zStak(:,kk),'linear','extrap');
            ageT(:,kk) = tmpAgeT;
        end
        pseudoAgeModel{ii} = ageT;
        pseudoAgeModel{ii}(1,:) = 0; % Set NaNs to Zero
        tmpModel = pseudoAgeModel{ii};
        clear('ageMod','zAxe','ageT');

        [m,n] = size(RadarStack{ii});
        q = 2; % Quantity times the current number of samples
        RadDeposit = zeros(q.*m,n);
        RadStack = RadarStack{ii};
        DepositionAxis{ii} = [linspace(min(tmpModel(:)),max(tmpModel(:)),length(tmpModel(:,1)))]';
        DepositAxe = DepositionAxis{ii};
        % ReSample in Age Domain
        DepositAxe = [linspace(min(DepositAxe),max(DepositAxe),q.*length(DepositAxe))]';
        DepositionAxis{ii} = DepositAxe;
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