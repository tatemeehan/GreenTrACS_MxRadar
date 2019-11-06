%% Time To Depth Conversion
    % Allocation
    RadarDepth = RadarStack;  
    
    for ii = 1:nFiles
        RadDepth = zeros(size(RadarStack{ii}));
        RadStack = RadarStack{ii};
        DepthAxe = DepthAxis{ii};
        zStak = zStack{ii};
        tStak = tStack{ii};
        parfor (kk = 1:size(RadarStack{ii},2), nWorkers)
            % [Time, Depth] - Array to be resampled
            tzcurve = [tStak(:,kk),zStak(:,kk)];
            % Perform Depth Conversion
            [RadDepth(:,kk)] = timeDepthConversion(RadStack(:,kk),tzcurve,DepthAxe);
        end
        RadarDepth{ii} = RadDepth;
    end
    clear('RadarStack','RadDepth','RadStack','DepthAxe','zStak','tStak');
