function [ W ] = nonParametricSmooth( x,y,xmod,R )
%nonParametricSmooth smooths a 1D data set using the biSquare Kernal 
%   Input:  x - independent variable
%           y - dependent variable
%           xmod - locations of estimates
%           R - Rank of Moving window
%
% Output: - W Smoothed data estimate
%
% Boise State University: Tate Meehan, GreenTrACS 2016
%% 
x = x(:);y=y(:);xmod = xmod(:);
W = zeros(size(xmod));
for ii = 1:length(xmod)
    dist = sqrt((x-xmod(ii)).^2);
    Ix = find(dist<=R);
    Ix = Ix(isfinite(y(Ix)));
    weight = 15/16*(1-(dist(Ix)./R).^2).^2;
    if isempty(weight)
        weight = NaN;
    end
    W(ii) = sum(weight(:).*y(Ix))./sum(weight);
end

end

