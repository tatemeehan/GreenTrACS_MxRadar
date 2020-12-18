function [gaindat]=gain(dat,tpow,rollOff)

%format [gaindat,tg]=gain(dat)

% Written by Bradford 
% Modified by Mikesell (01/19/2016)
% Meehan added rollOff 8/19 

nt = length( dat(1,:) );
ns = length( dat(:,1) );
if nargin < 3
    tg = [1:ns].^tpow.';
else
    if 2.*rollOff < ns
        rollOff = round(ns./2)+1;
        tg1 = [1:rollOff].^tpow.';
        tg2 = flipud(tg1);
        tg = [tg1;tg2(1:(ns-rollOff))];
    elseif ns > rollOff
        tg1 = [1:rollOff].^tpow.';
        tg2 = flipud(tg1);
        tg = [tg1;tg2(1:(ns-rollOff))];
    else
        tg = [1:ns].^tpow.';
    end
end
k=0;
gaindat=zeros(size(dat));
while (k < nt) % loop through each trace and apply
	gaindat(:,k+1)=dat(:,k+1).*tg;
k=k+1;
end
	
