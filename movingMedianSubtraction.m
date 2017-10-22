function [ dataOut ] = movingMedianSubtraction( data, R )
%   movingMedianSubtraction calculates the median of the nearest 2R+1
%   columns. A median Value is computed for each column. The Median Valued
%   matrix has dimensions of the data. The Median Matrix is subtraced from
%   the input matrix and is output as the result.
%   Input:  data - Matrix of Traces
%           R - Rank of Moving window
%
% Output: - dataOut median subtracted data
%
% Written by: Tate Meehan, Boise State University, GreenTrACS 2017

% Handle R Length Condition
if (2*R)+1 > size(data,2)
    R = floor((size(data,2)-1)./2);
end

x = [1:size(data,2)]; x = x(:); % Create Indexing Vector
M = zeros(size(data));          % Allocate Median Matrix
lag = round(2*R./10)+1;         % Overlapping Gates
lagIx = 1:lag:length(x);        % Gate Index
dist = zeros(length(lagIx),size(data,2));   % Allocate Distance Array

% Estimate Median within Gate for Each Lag Separation
for ii = 1:lag:length(x)
    dist(ii,:) = sqrt((x-x(ii)).^2);
    Ix = find(dist(ii,:)<=R);
    M(:,ii) = nanmedian(data(:,Ix),2); 
end

% Overlap Lagged Median Estimates
for ii = 1:lag:length(x)
    overlapDist = find(dist(ii,:)<=(lag/2)+1);
    M(:,overlapDist) = repmat(M(:,ii),1,length(overlapDist));
end
dataOut = data - M;

end