function [GPSix,edgeIx] = removeStaticGPS(Lon,Lat,win)
% removeStaticGPS algorithmically removes stationary GPS locations from ROV
% GPS data. The derivative of the variance of the GPS trace is classified
% into static and roving indecies by absolute mean deviation.

isDTM = 1; % Use Derivative Threshohold Method
isAMD = 0; % Use Absolute Deviation of Means Classification

% S the discretization for differentiaion must be an even valued integer.
win = round(win);
if mod(win,2) ~= 0
    S = win-1;
else
    S = win;
end

%% Derivative Threshold Method
% 800 times faster and Capable of Finding Start and Stop Indicies
if isDTM
%     tic
    % First Derivative
    sLon = abs((Lon((1*S+1):end)-Lon(1:(end-1*S)))./(1*S+1)).^1;
    sLat = abs((Lat((1*S+1):end)-Lat(1:(end-1*S)))./(1*S+1)).^1;
    % Append Edges to Maintain Indicies
    sLon = [sLon(1).*ones(floor((1*S+1)/2),1);sLon;sLon(end).*ones(floor((1*S+1)/2),1)];
    sLat = [sLat(1).*ones(floor((1*S+1)/2),1);sLat;sLat(end).*ones(floor((1*S+1)/2),1)];
    % Cut off Threshold (May not be Robust)
    threshold = 1e-6;
    GPSixLat = find(sLat > threshold);
    GPSixLon = find(sLon > threshold);
    
    % Take fewer Indicies Home
    if(length(GPSixLat)>length(GPSixLon))
        GPSix = GPSixLon(ismember(GPSixLon,GPSixLat));
    elseif(length(GPSixLat)<length(GPSixLon))
        GPSix = GPSixLat(ismember(GPSixLat,GPSixLon));
    else
        GPSix = GPSixLon(ismember(GPSixLon,GPSixLat));
    end
    % Take More Indicies Home
%     if(length(GPSixLat)<length(GPSixLon))
%         GPSix = GPSixLon(ismember(GPSixLon,GPSixLat));
%     elseif(length(GPSixLat)>length(GPSixLon))
%         GPSix = GPSixLat(ismember(GPSixLat,GPSixLon));
%     else
%         GPSix = GPSixLon(ismember(GPSixLon,GPSixLat));
%     end
    
    % First Derivative of Good Indicies
    sGPSix = abs((GPSix((1*S+1):end)-GPSix(1:(end-1*S)))./(1*S+1)).^1;
    % Append Edges to Maintain Indicies
    sGPSix = [sGPSix(1).*ones(floor((1*S+1)/2),1);sGPSix;sGPSix(end).*ones(floor((1*S+1)/2),1)];
    % Slopes Greater than 1 Will be Static Points
    s1GPSix = find(sGPSix > 1);
    % Difference Static Indicies and Append with 1 to keep Index Correct
    ds1GPSix = [1;diff(s1GPSix)];
    % Bin the Data (Should be N .nc Files)
    bins = [1;sort([find(ds1GPSix > .1.*max(ds1GPSix));find(ds1GPSix > .1.*max(ds1GPSix))-1]);length(ds1GPSix)];
    % Loop Through Bins to find Midpoint, More Accurate bc Slope Smear
    for ii = 1:2:length(bins)
        avgIx(ii) = round(mean(s1GPSix(bins(ii:ii+1))));
    end
    avgIx(avgIx == 0) = [];
    % The Start and Stop Indicies of the GPS/GPR Files.
    edgeIx = sort([avgIx(:);avgIx(:)-1]);
    edgeIx(1) = [];
    edgeIx = [edgeIx;length(GPSix)];
    edgeIx = reshape(edgeIx,[2,length(edgeIx)/2]);
%     toc
end
%% Absolute Deviation from Mean Method
% Good at Removing Static Points, Unable to Seek Breaks in GPR Files
if isAMD
    % Compute Running Variance of Lon and Lat for Window Size N
    varLat = movvar(Lat,win);
    varLon = movvar(Lon,win);
    
    
    % Compute Qunatile Function
    [qvarLat, latIx] = sort(varLat);
    [qvarLon, lonIx] = sort(varLon);
    % Compute Running Differenced Average
    qvarLatAvg = (qvarLat((1*S+1):end)-qvarLat(1:(end-1*S)))./(1*S+1);
    % Append Edges to Maintain Indicies
    qvarLatAvg = [qvarLatAvg(1).*ones(floor((1*S+1)/2),1);qvarLatAvg;qvarLatAvg(end).*ones(floor((1*S+1)/2),1)];
    % Transform Slope
    qvarLatS = qvarLatAvg.^.1;
    % Normalize Slope
    % qvarLatS = qvarLatAvg./max(qvarLatAvg);
    
    qvarLonAvg = (qvarLon((1*S+1):end)-qvarLon(1:(end-1*S)))./(1*S+1);
    % Append Edges to Maintain Indicies
    qvarLonAvg = [qvarLonAvg(1).*ones(floor((1*S+1)/2),1);qvarLonAvg;qvarLonAvg(end).*ones(floor((1*S+1)/2),1)];
    % Transform Slope
    qvarLonS = qvarLonAvg.^.1;
    % Normalize Slope
    % qvarLonS = qvarLonAvg./max(qvarLonAvg);
    
    % Optimize Mean Difference of Classifications
    qLon = 0.95;
    % qLat = 1;
    kk = 1;
    atol = 30;
    absErr = atol+1;
    rmIx = 0;
    tic
    while absErr > atol || abs(rmIx - length(qvarLonS)) < 0.05.*length(qvarLonS)
        for ii = round(length(qvarLonS) - qLon.*length(qvarLonS))+1:round(qLon.*length(qvarLonS)-1)
            tmpLon = abs(mean(qvarLonS(1:ii)) - mean(qvarLonS(ii+1:end)));
            tmpLat = abs(mean(qvarLatS(1:ii)) - mean(qvarLatS(ii+1:end)));
            if ii == round(length(qvarLonS) - qLon.*length(qvarLonS))+1
                rmIxLon(kk) = ii;
                testLon = tmpLon;
            elseif tmpLon > testLon
                rmIxLon(kk) = ii;
                testLon = tmpLon;
            end
            if ii == round(length(qvarLonS) - qLon.*length(qvarLonS))+1
                rmIxLat(kk) = ii;
                testLat = tmpLat;
            elseif tmpLat > testLat
                rmIxLat(kk) = ii;
                testLat = tmpLat;
            end
        end
        
        absErr = abs(rmIxLat(kk)-rmIxLon(kk));
        testErr(kk) = absErr;
        % rmIx = round(mean([rmIxLat(kk),rmIxLon(kk)]));
        rmIx = round(min([rmIxLat(kk),rmIxLon(kk)]));
        
        if absErr > 30 || abs(rmIx - length(qvarLonS)) < 0.05.*length(qvarLonS)
            kk = kk +1;
            qLon = qLon - 0.025;
            %     qLat = qLat- 0.025;
            if qLon < 0.9
                [~,minIx] = min(testErr);
                kk = minIx;
                break
            end
            
        end
        
    end
    toc
    
    
    lonIx(1:rmIx) = [];
    % latIx(1:rmIx) = [];
    
    GPSix = sort(lonIx);
end

end

