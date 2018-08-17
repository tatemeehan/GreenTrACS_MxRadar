    %% Partial Stack
    % Selected Offsets 4 and 12 m for Stack. The Late Times of 4m are muted
    RadarStack = cell(nFiles,1);
    [~,ix4] = min(abs(offsetArray - 4));
    [~,ix12] = min(abs(offsetArray - 12));
    for ii = 1:nFiles
        % Bottom Mute
        stackTaper = zeros(size(RadarNMO{ix4,ii},1),1);
        stackTaper(1:200,1) = tukeywin(size(RadarNMO{ix4,ii}(1:200,:),1),0.05);
        stackTaper = stackTaper*ones(1,size(RadarNMO{ix4,ii},2));
        % Stack Offsets
        RadarStack{ii} = (RadarNMO{ix4,ii}.*stackTaper)+RadarNMO{ix12,ii};
        clear stackTaper
    end