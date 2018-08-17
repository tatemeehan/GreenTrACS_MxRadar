function [ V ] = movingVariance( data, R, dim, w )
%   movingVariance calculates the variance of the nearest 2R+1
%   neighbors of data dimension dim. data can be 1 or 2-D. 
%   Input:  data - Vector or Matrix of spatial - temporal data
%           R - Rank of Moving window
%           dim - dimension to compute variance
%           w - weighting: 0 is sample variance 1 is population variance
%                          w may be a vector of weights length data in the
%                          chosen dimension.
%
% Output: - V data variance
%
% Written by: Tate Meehan, Boise State University, GreenTrACS 2017

% Handle Data Conditions
[nr, nc] = size(data);

% If only data is supplied
if nargin < 2
    % Default Window is 5 Points
    R = 2;
    % Default Weigting is None
    w = 0;
    % Search for Dimensionality 
    if nr == 1
        dim = 2;
    elseif nc == 1
        dim = 1;
    else
        dim = 2;
    end
end
% If data and R is supplied
if nargin < 3
    % Default Weigting is None
    w = 0;
    % Search for Dimensionality    
    if nr == 1
        dim = 2;
    elseif nc == 1
        dim = 1;
        
    else
        dim = 2;
    end
end
% If data, R, and dim is supplied
if nargin < 4
    % Default Weigting is None    
    w = 0;
end
%% Allocate Memory
x = [1:size(data,dim)]; x = x(:);       % Create Indexing Vector
if dim == 1
    V = zeros(1,size(data,dim));        % Allocate Median Matrix
else
    V = zeros(size(data,dim),1);
end
%% Moving Variance Computation
for ii = 1:length(x)
    dist = sqrt((x-x(ii)).^2);
    Ix = find(dist<=R);
    if dim == 1
        V(ii) = nanvar(data(Ix,:),w,dim);
    else
        V(ii) = nanvar(data(:,Ix),w,dim);
    end
end

end