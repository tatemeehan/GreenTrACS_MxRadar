function [ removeTraces, staticTraces ] = removeStaticTrace( Rad, multiplexNtrcs, nearChan,  nChan )
% removeCovTrace seeks and removes the traces that are duplicated when
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

% Surpress Warning 
warning('off','signal:findpeaks:largeMinPeakHeight')

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

% Deteminie Quantile Threshold by Maximum Difference of Mean Covariance
quantiles = linspace(0.985,1,1000);
meanMovingCov = zeros(1000,1);meanStaticCov = meanMovingCov;
for kk = 1:1000

qCov = quantile(qCovAvg(100:end),quantiles(kk));
tmp(kk) = qCov;

% Find Peak Change in Covariance - Evaluated q Probability
[~,peakIx] = findpeaks(qCovAvg(100:end),'MinPeakHeight',qCov);
dumTraces = traces(100:end);
peakIx = dumTraces(peakIx);

% Evaluate the Data Covariance Threshold
[Threshold] = min(peakIx); 

meanMovingCov(kk) = mean(covRad(Ix(1:Threshold-1)));
meanStaticCov(kk) = mean(covRad(Ix(Threshold:end)));

end

% Maximization Classification of Quantile Means
% [~,qIx] = max(diff(meanMovingCov));
% [~,qIx] = max(meanMovingCov(11:end)-meanMovingCov(1:end-10));

[~,qIx] = max(abs(meanStaticCov-meanMovingCov));


% q is Optimal Quantile
q = quantiles(qIx);

% Define Peak Covarince Threshold of Running Differenced Average Covariance
% Steer Clear of Padded Edge (100:end) - Creates Artificial Peak
qCov = quantile(qCovAvg(100:end),q);

% Find Peak Change in Covariance - Evaluated at Optimal Probability
[~,peakIx] = findpeaks(qCovAvg(100:end),'MinPeakHeight',qCov);
dumTraces = traces(100:end);
peakIx = dumTraces(peakIx);

% Evaluate the Data Covariance Threshold
[Threshold] = min(peakIx);

% Identify Static Traces
staticTraces = Ix(Threshold:end);

% Refinement Step
% We sometimes wish to add a few more static points or remove a few
% Simple Method of Adding 20 STD to mean and
meanMovingCov = mean(covRad(Ix(1:Threshold-1)));
stdMovingCov = std(covRad(Ix(1:Threshold-1)));

meanStaticCov = mean(covRad(Ix(Threshold:end)));
% Weighted average of moving and static covariances
ThresholdCov = mean([meanMovingCov,meanStaticCov,meanMovingCov,meanMovingCov,meanMovingCov]);
staticTraces = find(covRad>ThresholdCov);

% Redux of Slope Thresholding
% % Compute Qunatile Function
% [qCovRad, ~] = sort(covRad(staticTraces));
% 
% % Compute Running Differenced Average
% qCovAvg = (qCovRad((1*S+1):end)-qCovRad(1:(end-1*S)))./(1*S+1);
% % Append Edges to Maintain Indicies
% qCovAvg = [qCovAvg(1).*ones(floor((1*S+1)/2),1);qCovAvg;qCovAvg(end).*ones(floor((1*S+1)/2),1)];
% 
% % Define Peak Covarince Threshold of Running Differenced Average Covariance
% % Steer Clear of Padded Edge (100:end) - Creates Artificial Peak
% qCov = quantile(qCovAvg(100:end),q);
% 
% % Find Peak Change in Covariance - Evaluated at Optimal Probability
% [~,peakIx] = findpeaks(qCovAvg(100:end),'MinPeakHeight',qCov);
% dumTraces = traces(100:end);
% peakIx = dumTraces(peakIx);
% 
% % Evaluate the Data Covariance Threshold
% [Threshold] = min(peakIx);
% 
% % [~,~,goodIx] = ismember(covRad(staticTraces),qCovRad(1:Threshold-1));
% staticTraces(1:Threshold-1) = [];

% % Order of Magnitude Filter in case some points sneak through
% % Remove Points
% tmp = meanStaticCov(qIx); mag = floor(log10(tmp));
% rmvIx = covRad(staticTraces)<1.*10.^mag;
% staticTraces(rmvIx) = [];
% goodTraces = traces;
% goodTraces(staticTraces) = [];
% % Add Points
% addIx = (goodTraces(covRad(goodTraces)>1.*10.^mag))';
% staticTraces = [staticTraces;addIx];

% Identify Multiplexed Static Traces
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