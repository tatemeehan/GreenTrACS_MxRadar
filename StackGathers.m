    %% Partial Stack
    % Selected Offsets 4, 9, and 12 m for Stack. 
    % The Late Times of 4m are muted
    RadarStack = cell(nFiles,1);
    [~,ix4] = min(abs(offsetArray - 4));
    [~,ix9] = min(abs(offsetArray - 9.33));
    [~,ix10] = min(abs(offsetArray - 10.67));
    [~,ix12] = min(abs(offsetArray - 12));
    taperWin = 1:200; % Time Samples for Taper Window
    for ii = 1:nFiles
        % Bottom Mute
        stackTaper = zeros(size(RadarNMO{ix4,ii},1),1);
        stackTaper(taperWin,1) = tukeywin(size(RadarNMO{ix4,ii}(taperWin,:),1),0.05);
        stackTaper = stackTaper*ones(1,size(RadarNMO{ix4,ii},2));
        % Stack Offsets
        RadarStack{ii} = (RadarNMO{ix4,ii}.*stackTaper)+RadarNMO{ix12,ii}...
            +RadarNMO{ix9,ii} + RadarNMO{ix10,ii};
        clear stackTaper
    end