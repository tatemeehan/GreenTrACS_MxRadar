function [ Stak ] = StakR( data, R, type )
% StakR stacks adjacent traces determined by R; the stacking window rank,
%   and type, a string identifying which  filter to be applied. The User
%   may apply sum, mean, or weighted mean stacking. 'kernel' the string
%   identifying a hamming weighted mean is applied by default.
%   Window Size is length 2R+1 and is centered about column(:,ii)
%
%   Inputs:     data            - Data Matrix
%               R               - Rank of Stacking Window
%               type            - 'sum' | 'mean' | 'kernel'  Stack Type
%
%   Outputs:    Stak            - Stacked Data Matrix
%
% Written by Tate Meehan, Boise State University, GreenTrACS 2016-2017

if nargin == 2      % Version Control, Default is Kernel Weighted Stacking
    type = 'kernel';
end

[nr, nc] = size(data);

% Appropiate Window Size
if R > nc/2+1
    R = floor(nc/2);
end
% Allocate Stacked Data Matrix
Stak = zeros(nr,nc);

% Stack Traces
% Sum Stack
if strcmp(type,'sum')
    for ii = 1: nc
        if ii <= R
            % Leading Stak Edge Taper
            Stak(:,ii) = (sum(data(:,1:ii+R),2))... Amplitude Normalize
                +((R-ii)./R).*((sum(data(:,1:ii+R),2))./nr);
        elseif ii > R && ii <= (nc-R)
            % Sum Central Traces
            Stak(:,ii) = sum(data(:,ii-R:ii+R),2);
        else
            % Trailing Stak edge Taper
            Stak(:,ii) = sum(data(:,ii-R:nc),2)... Amplitude Normalize
                +(1-(length(data(:,ii-R:nc))./(2*R))).*(sum(data(:,ii-R:nc),2)./nr);
        end
    end
end

%  Mean Stack
if strcmp(type,'mean')
    for ii = 1: nc
        if ii <= R
            % Leading Stak Edge Taper
            Stak(:,ii) = (mean(data(:,1:ii+R),2));
        elseif ii > R && ii <= (nc-R)
            % Mean Central Traces
            Stak(:,ii) = mean(data(:,ii-R:ii+R),2);
        else
            % Trailing Stak edge Taper
            Stak(:,ii) = mean(data(:,ii-R:nc),2);
        end
    end
end

% Weighted Mean Stack
if strcmp(type,'kernel')
    kernel = hamming(2.*R +1);
    Stak = conv2(1,kernel./sum(kernel),data(:,:),'same');
    
end
end

    