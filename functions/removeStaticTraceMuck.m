function [ removeTraces, staticTraces ] = removeStaticTraceMuck( Rad, multiplexNtrcs, nearChan,  nChan )
% removeCovTrace seeks and removes the traces which are duplicated when
% the radar array is stationary. This only applies to FreeRun DAQ. 
% This approach uses the covariance matrix to self determine an appropriate
% threshold for static trace removal. 
% Adapted for Single Channel Data.
%
%   Inputs  Rad   - The Near Offset Channel Data (Single Channel Data)
%           ~~ For MultiChannel Data  ~~
%           nChan - The Number is Radar Channels
%           plexTraces - The Multiplexed Trace Indicies
%           nearChan - The Channel Number of the Near Offset Data 
%
%   Output removeTraces - Static Trace indicies
%
% Written by, Tate Meehan, Boise State University, GreenTrACS 2017

[~, ntrc] = size(Rad); % Evaluate Size of Data
traces = 1:ntrc; % Create Array of Trace Indicies

% For Single Channel Data
if nargin == 1
    nChan = 1;
    nearChan = 1;
    multiplexNtrcs = ntrc;
end

modChan = nChan - nearChan; % Final Multiplexed Channel Indicies
nearPlex = nearChan:nChan:multiplexNtrcs; % Multiplexed Near-Channel Indicies

% Compute Covariance of Traces Distance S Apart
S = 10; 
covRad = diag(cov(Rad),S);

% Infill Truncated Diagonal Elements
covRad = [ones(S/2,1).*covRad(1);covRad;ones(S/2,1).*covRad(end)];

% Compute Qunatile Function
[qCovRad, Ix] = sort(covRad);

% Compute Running Differenced Average
qCovAvg = (qCovRad((1*S+1):end)-qCovRad(1:(end-1*S)))./(1*S+1);
% Append Edges to Maintain Indicies
qCovAvg = [qCovAvg(1).*ones(floor((1*S+1)/2),1);qCovAvg;qCovAvg(end).*ones(floor((1*S+1)/2),1)];

% Compute 99% Quantile of Running Differenced Average to Define Peak Thresh
% Steer Clear of Padded Edge (100:end) - Creates Artificial Peak
q99 = quantile(qCovAvg(100:end),0.99);

% Find Peak Change in Covariance - Evaluated 99% Probability
% [~,peakIx] = findpeaks(qCovAvg(100:end),traces(100:end),'MinPeakHeight',q99);
[~,peakIx] = findpeaks(qCovAvg(100:end),'MinPeakHeight',q99);
dumTraces = traces(100:end);
peakIx = dumTraces(peakIx);

% Evaluate the Data Covariance Threshold
[Threshold] = min(peakIx); threshCov = qCovRad(Threshold);

% Identify Static Traces
staticTraces = Ix(Threshold:end);

% Sort Static Traces and Find Distances
sortStaticTraces = sort(staticTraces);
group = diff(sortStaticTraces);

% Remove Traces During Array Turning
groupIx = find(group > 1 & group < 250);

% Append Additional Traces for Removal
for jj = 1:length(groupIx)
    appendTraces = sortStaticTraces(groupIx(jj)) + 1 : sortStaticTraces(groupIx(jj) + 1) -1;
    staticTraces = [staticTraces;appendTraces(:)];
end

% staticTraces = find(covRadW > Threshold ); % Static Traces
staticPlexTraces = nearPlex(staticTraces) + modChan; 

% Infill Multiplex Trace Indicies For Removal
removeBin = zeros(length(staticPlexTraces),nChan);
for ii = 1:nChan
    if ii == nChan
        removeBin(:,ii) = staticPlexTraces;
    else
    removeBin(:,ii) = staticPlexTraces - (nChan - ii);
    end
    % Indicies of Static Traces
    removeTraces = removeBin(removeBin(:) <= multiplexNtrcs);
end

end

