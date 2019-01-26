%% Normal Moveout Correction
    % Maximum Allowable Wavelet Stretch [%]
    Stretch = 100;
    % Allocation
    xStack = cell(nFiles,1);
    tStack = cell(nFiles,1);
    vStack = cell(nFiles,1);
    zStack = cell(nFiles,1);
    DepthAxis = cell(nFiles,1);
    DepthMatrix = cell(nFiles,1);
    RadarNMO = Radar;
    for ii = 1:nFiles
        jj = 1;
        isInterpolate = 1;
        % Perform NMO Correction with Grid Interpolation
        [RadarNMO{jj,ii},xStack{ii},tStack{ii},vStack{ii},~] = ...
            commonOffsetNMO(Radar{jj,ii},dt,f0,offsetArray(jj),TraverseX{ii},StackingTime{ii},StackingVelocity{ii},Stretch,isInterpolate);
        % Grid of Depths
        zStack{ii} = (vStack{ii}.*tStack{ii})./2;
        % Axis for Depth Image
        DepthAxis{ii} = (linspace(min(zStack{ii}(:)),max(zStack{ii}(:)),size(RadarNMO{jj,ii},1)))';
        DepthMatrix{ii} = DepthAxis{ii}(:)*ones(1,size(RadarNMO{jj,ii},2));
        
        % Allocation for Parfor Overhead
        dumX = xStack{ii};
        dumT = tStack{ii};
        dumV = vStack{ii};
        parfor (jj = 2:nChan, nWorkers)
%         for jj = 2:nChan
            % Perform NMO Correction on Remaining Channels w/o Interpolation
            isInterpolate = 0;
            [RadarNMO{jj,ii},~,~,~,~] = ...
                commonOffsetNMO(Radar{jj,ii},dt,f0,offsetArray(jj),dumX,dumT,dumV,Stretch,isInterpolate);
        end
        clear('dumX','dumT','dumV');
    end