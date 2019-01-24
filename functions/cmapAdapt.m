function [c] = cmapAdapt(data,cmap)
% cmapAdapt shifts the zero crossing of the divergent colormap to the
% median of the data. 
%   The user specifies the data and their choice of colormap.
%   MATLAB defalt colormaps may be specified as a string or otherwise.
%
% Boise State University: Tate Meehan, NASA ISGC 2019

% Test if input is a MATLAB default string
if ischar(cmap)
    cmap = colormap(cmap);
end
n = length(cmap);
%Normalize Data
data = 2.*((data(:)-min(data(:)))./(max(data(:))-min(data(:))))-1;
% Compute Min,Median, and Max of Data
dq = [quantile(data(:),[0,.5,1]),floor(quantile(1:n,.5))];
rangedq = dq(3)-dq(1);
lhs = rangedq./n;
% Shift Median to Zero
rhs = dq(2);
% Compute number of Indicies to Shift Color Map
ix = round(rhs./lhs);
% Shift Colormap Center and Interpolate 
a = [interp1(linspace(1,dq(4)+ix,dq(4)),cmap(1:dq(4),:),1:(dq(4)+ix))];
b = [interp1(linspace(dq(4)+ix+1,n,n-dq(4)),cmap(dq(4)+1:end,:),(dq(4)+ix+1):n)];
c = [a;b];
end

