% Multiple S&S GPR Files are Concatenated into a single radar image.
% Cross Correlation is Performed to Optimize the Alignment of the Images

    % Grab Early Time Data For AirWave Cross Correlation Alignment
    % 2016 Approximate Direct Wave Arrival Times
    if strcmp(Year{1},'2016')
        xcorrWindow = [30, 105, 180, 1, 80, 145, 1, 45, 130;...
                    60, 150, 210, 40, 110, 175, 5, 80, 160] + padding;        
%         xcorrWindow = [30, 105, 180, 0, 80, 145, -20, 45, 130;...
%                     60, 150, 210, 40, 110, 175, 5, 80, 160] + padding;
    end
    % 2017 Approximate Direct Wave Arrival Times
    if strcmp(Year{1},'2017')
%         xcorrWindow = [175, 250, 85, 140, 230, 75, 100, 180, 25;...
%                     200, 300, 150, 160, 260, 115, 150, 220, 75] + padding;
        xcorrWindow = [90, 155, 220, 65, 125, 195, 45, 110, 175;...
                    115, 190, 260, 90, 165, 265, 75, 135, 205] + padding;
    end
    if isKill
        xcorrWindow(:,killChan) = [];
    end
    
    % Cross-Correlation for multipleFile Alignment Prior to Concatenation
    Radar = xcorrAlignTraces( Radar, xcorrWindow );
    
    % Concatenate Files
    isCat = 1;
    if isCat
        nFiles = 1;
        TraverseDistance = sum(TraverseDistance);
        catRadar = cell(nChan,nFiles);
        for jj = chan            
            catRadar{jj} = [Radar{jj,:}];
        end
         Radar = catRadar; clear catRadar;
         tmptrhd = [trhd{:}];
         trhd = cell(1,nFiles);
         trhd{nFiles} = tmptrhd;
         TimeAxis(nFiles+1:end) = [];
         tmpXYZ = cell2mat(XYZ);
         XYZ = tmpXYZ;
         clear('tmptrhd','tmpXYZ');
    end