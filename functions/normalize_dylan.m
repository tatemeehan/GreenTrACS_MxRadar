function a = normalize_dylan(a)
%
% Normalize traces by max(abs(trace)).
%
% USAGE: a = normalize_dylan(a)
%
% INPUT:
%   a(npts,ntrc)=data matrix
%
% OUTPUT:
%   a(npts,ntrc)=trace normalized data matrix
%
% Created by Dylan Mikesell: 4 March 2013

[~,ntrc]=size(a);

for ii=1:ntrc
    % check to make sure this is not a zero trace (we don't want to touch null traces)
    if( sum(a(:,ii))~=0 )
        a(:,ii)=a(:,ii)./max(abs(a(:,ii))); % normalize by largest amplitude in trace
    end
    
end

return