function [Rad] = trcNormalize(Rad)
% Trace Normalize reduces the median of the radargram to zero and ensures
% that its amplitudes are between -1 and 1
tmp = 2.*((Rad(:)-min(Rad(:)))./(max(Rad(:))-min(Rad(:))))-1;
dq = quantile(tmp(:),[0.425,.5,.575]);
tmp = tmp - dq(2);
% reMap Min and Max values to -1 and 1
maxIx = find(tmp>1);
minIx = find(tmp<-1);
if isempty(maxIx)
    [~, maxIx] = max(tmp);
    tmp(maxIx) = 1;
else
    tmp(maxIx) = 1;
end

if isempty(minIx)
    [~, minIx] = min(tmp);
    tmp(minIx) = -1;
else
    tmp(minIx) = -1;
end
% RadarStack is Normalized and Reduced to Zero Median
Rad = reshape(tmp,size(Rad));
clear('tmp')
end

