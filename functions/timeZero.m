function [ timeZeroData, twt, to ] = timeZero( data, t, dt, R, pow, offset, isCMP )
% timeZero identifies the first arrival time and trims negative time from 
%   the data.  This function is used as the "time zero" correction. And
%   assumes time zero to be identical across all traces.
%
%   Input:      data - Uncorrected Radargram Data Matrix
%               t    - Two-way Travel Time Window   
%               R    - Rank of the MER Window [ns]
%               pow  - Power Function of MER Picker 0, 1, 2, 3
%   Output:     timeZeroData - Time Zero Correced Data
%               timeZero     - Zero Corrected TwoWayTime
%               to           - Intercept Time (Removed Time)
%
% Written by Tate Meehan, Boise State University, GreenTrACS 2016

% Error Handeling
if nargin < 7
    isCMP = 0;
elseif nargin < 6
    isCMP = 0;
end
% Allocation for MER Picks
CommonFBPick = zeros(size(data,2),1);
% Modified Energy Ratio First Break Picking
if isCMP
    % Find Minimum Offset Index
    [~, nearOffsetIx] = min(offset);
    % Set Minimum Offset
    offset = offset(nearOffsetIx);
    % Grab Near offset Trace
    tmp = data(:,nearOffsetIx);
    % First Break Picking
    [CommonFBix, ~, ~] = wong_mer(tmp, R, dt, pow);
    CommonFBPick = t(CommonFBix);
else
    for kk = 1:length(data(1,:))
        % First Break Picking
        tmp = data(:,kk);
        [CommonFBix, ~, ~] = wong_mer(tmp, R, dt, pow);
        CommonFBPick(kk) = t(CommonFBix);
    end
end
% Travel Time Correction is the Mean Static Adjustment
FBpick = mean(CommonFBPick);
% Find Good Travel Times            
timeIx = find((t >= abs(FBpick)));
% Non-Zero Instrument Noise For Padding
pad = mean(data(1:timeIx-1,:));
% Truncate Data at Airwave Arrival
timeZeroData = data(timeIx,:);

% Pad Data with Instrument Noise for Non-zero offset Data
c = 0.3; %[m/ns]
padsize = round(offset./c./dt);
timeZeroData = padarray(timeZeroData,padsize,mean(pad),'pre');

% Find Good Travel Times and Time Zero from MER Picks
tmpTime = t(timeIx);
to = tmpTime(1);
% Pad Travel Time axis
padto = to - (padsize:-1:1).*dt;
tmpTime = [padto(:);tmpTime];
% Re-Calculate Time Zero
to = tmpTime(1);
% Re-Configure Travel Time Array
twt = tmpTime-to;

end

