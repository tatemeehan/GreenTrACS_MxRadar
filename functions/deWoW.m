function [filtdat]=deWoW(dat,dt,f0)
tic
%written by John Bradford 6/22/2005
%modified by Bradford from Liberty's code 6/24/2005
% Modified by Tate Meehan, Adjusted Filter Length lf, GreenTrACS 2016
%applies low cut filter to data - uses sensors and software moving boxcar 
%filter
%format dewow(dat,dt,f0)


nt=length(dat(1,:));
ns=length(dat(:,1));
% Convert Units
f0Hz = f0 * 1e6;       % [Hz]
dtSec = dt * 1e-9;      % [s]
lf=ceil(1/f0Hz/dtSec);
bg=zeros(ns-lf-1,nt);
filtdat=zeros(size(dat));
bgs=mean(dat(1:2*lf+1,:));
bge=mean(dat(ns-2*lf:ns,:));
for p=1:ns
	if(p < lf+1 )
		filtdat(p,:)=dat(p,:)-bgs;
	elseif (p > ns-lf-1)
		filtdat(p,:)=dat(p,:)-bge;
	else
		bg(p,:)=mean(dat(p-lf:p+lf,:));
		filtdat(p,:)=dat(p,:)-bg(p,:);
	end
end
toc
		