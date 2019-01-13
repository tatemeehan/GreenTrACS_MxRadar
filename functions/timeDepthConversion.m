function [ztrc]=timeDepthConversion(trc,tzcurve,z)
% [ztrc]=timeDepthConversion(trc,tzcurve,z)
% Original by G.F. Margrave (1995)
% [ztrc,z]=time2depth(trc,t,tzcurve,dz)
% TIME2DEPTH: Convert a trace from time to depth by a simple stretch
%
% [ztrc,z]=time2depth(trc,t,tzcurve,dz)
%
% TIME2DEPTH converts a single trace from time to depth. The conversion is
% specified by a time-depth curve that is stored in a two-column matrix.
% The first column is a list of times and the second is a list of 
% corresponding depths. The times can be either one-way or two-way and
% need only be consistent with the time-coordinate vector of the input.
% Between points in the time-depth curve, depths are interpolated linearly.
% Note that this function does not apply an antialias filter. It is
% up to the user to ensure that dz is sufficiently small to preclude
% aliasing.
% 
% trc ... the trace in time
% t ... time coordinate vector for trc
% tzcurve ... an n-by-2 matrix giving the time-depth curve.
%	n is the number of points on the curve and is arbitrary.
%	(more points is usually more accurate.) The first column is
%	time and the second column is the depths that correspond to
%	the times in the first column. The first time should be zero
%	and the last should be greater than or equal to max(t). This can be
%	created with sonic2tz.
% dz ... depth sample size. 
%
% NOTE: to avoid aliasing pick dz<vmin*dt/n where n is 1 for one-way time
% and 2 for 2way time and vmin is the slowest velocity in the model.
%
% G.F. Margrave, CREWES, Nov, 2000
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its 'AUTHOR' (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) This SOFTWARE may be used by any individual or corporation for any purpose
%    with the exception of re-selling or re-distributing the SOFTWARE.
%
% 2) The AUTHOR and CREWES must be acknowledged in any resulting publications or
%    presentations
%
% 3) This SOFTWARE is provided "as is" with no warranty of any kind
%    either expressed or implied. CREWES makes no warranties or representation
%    as to its accuracy, completeness, or fitness for any purpose. CREWES
%    is under no obligation to provide support of any kind for this SOFTWARE.
%
% 4) CREWES periodically adds, changes, improves or updates this SOFTWARE without
%    notice. New versions will be made available at www.crewes.org .
%
% 5) Use this SOFTWARE at your own risk.
%
% END TERMS OF USE LICENSE

% time2depth Modified by: Tate Meehan, Boise State University, GreenTrACS
% Rather than seeking z1 and z2 and providing dz for the creating depth 
% axis for trace interpolation, the deisred depth axis is input.
% depthTimeConversion.m is a wrapper for all subroutines of Margrave's
% time2depth.m

% Define Time Axis
t = tzcurve(:,1);

%input trace must be regularly sampled
if(sum(abs(diff(diff(t))))>.00001)
    error('input time trace must be sampled regularly');
end


%check tz curve limits
tz=tzcurve(:,1);zt=tzcurve(:,2);
if(min(tz)~=0)
	error('tz curve must start at zero time');
end
if(max(tz)<max(t))
	error('tz curve must extend to times greater than max(t)');
end

%make sure depths are monotonic
ztest=diff(zt);
ind=find(ztest<=0, 1);
if(~isempty(ind))
    zt = sort(zt);
% 	error('depths on tzcurve must increase monotonically')
end

%interpolation sites
tint=pwlint(zt,tz,z);

%sinc function interpolation
ztrc=sinci(trc,t,tint);
end
function yi=pwlint(x,y,xi)
% PWLINT: piecewise linear interpolation (much faster than interp1)
%
% yi=pwlint(x,y,xi)
%
% PWLINT performs linear interpolation when it is known that that
% input function (x,y) is piecewise linear. If length(x) is much less than
% the length(xi), this is MUCH faster than the built in INTERP1. Points
% in xi which are outside the bounds of x will return nans.
%
% G.F. Margrave
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its 'AUTHOR' (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) This SOFTWARE may be used by any individual or corporation for any purpose
%    with the exception of re-selling or re-distributing the SOFTWARE.
%
% 2) The AUTHOR and CREWES must be acknowledged in any resulting publications or
%    presentations
%
% 3) This SOFTWARE is provided "as is" with no warranty of any kind
%    either expressed or implied. CREWES makes no warranties or representation
%    as to its accuracy, completeness, or fitness for any purpose. CREWES
%    is under no obligation to provide support of any kind for this SOFTWARE.
%
% 4) CREWES periodically adds, changes, improves or updates this SOFTWARE without
%    notice. New versions will be made available at www.crewes.org .
%
% 5) Use this SOFTWARE at your own risk.
%
% END TERMS OF USE LICENSE

if(length(x)<=length(xi))
    nsegs=length(x)-1;
    yi=nan*zeros(size(xi));

    for k=1:nsegs

        %find the points in this line segment
        ii=between(x(k),x(k+1),xi,2);

        if( ii )
            % interpolate
            yi(ii)=y(k)+(y(k+1)-y(k))*(xi(ii)-x(k))/(x(k+1)-x(k));
        end

    end
else
    yi=nan*zeros(size(xi));
    for k=1:length(xi)
        ii=surround(x,xi(k));
        if(~isempty(ii))
            yi(k) = y(ii)*(x(ii+1)-xi(k))/(x(ii+1)-x(ii))+y(ii+1)*(xi(k)-x(ii))/(x(ii+1)-x(ii));
        else
            yi(k)=nan;
        end
    end
end
end
function trout=sinci(trin,t,tout,sizetable)
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its 'AUTHOR' (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) This SOFTWARE may be used by any individual or corporation for any purpose
%    with the exception of re-selling or re-distributing the SOFTWARE.
%
% 2) The AUTHOR and CREWES must be acknowledged in any resulting publications or
%    presentations
%
% 3) This SOFTWARE is provided "as is" with no warranty of any kind
%    either expressed or implied. CREWES makes no warranties or representation
%    as to its accuracy, completeness, or fitness for any purpose. CREWES
%    is under no obligation to provide support of any kind for this SOFTWARE.
%
% 4) CREWES periodically adds, changes, improves or updates this SOFTWARE without
%    notice. New versions will be made available at www.crewes.org .
%
% 5) Use this SOFTWARE at your own risk.
%
% END TERMS OF USE LICENSE

% SINCI is a wrapper for interpbl

if(nargin==4)
    n=sizetable(1);
    trout=interpbl(t,trin,tout,round(n/2));
else
    trout=interpbl(t,trin,tout);
end
end
function sint=interpbl(t,s,tint,n,m)
%INTERPBL ... bandlimited sinc function interpolation
%
% sint=interpbl2(t,s,tint,n,m)
%
% The sampling theorem gives the reconstruction formula to recover the
% continuous signal from its samples provided that the continuous signal is
% bandlimited (i.e. has compact support in the frequency domain). This
% formula is sint(k)=sum(sinc.*s) where sinc is a properly positioned sinc
% function (it must be shifted to have its maximum at the time of the
% interpoolated sample), s are the samples, and sint is the interpolated
% value. The sum is theoretically over an infinite set of samples and the
% sample weights from from the sinc function. To make this practical, the
% sinc function must be truncated so that the sum becomes finite. The sinc
% should also be smoothly tapered. Here, the half-length of the sinc
% function is specified by the input argument 'n' meaning the sinc extends
% for n points on either side of the maximum. So, if n=4, we say we are
% doing an 8 point interpolation. The truncated sinc is also windowed by a
% Gaussian whose standard deviation is n/2 so that the Gaussian is two
% standard deviations down at the truncation points.
%
% t ... vector of times of the input samples. t must be regularly spaced.
% s ... input samples from the bandlimited signal(s). This can be one trace
%       or a gather of traces. If the latter, then the interpolation sites
%       must be the same on all traces in the gather.
% *** the length of t must match one dimension of s. Preferably the row
% *** dimensions of both are the same.
% tint ... vector of times at which the interpolated values are desired
%       (the interpolation sites). tint does not need to be regular
%       although it can be.
% n ... half-length of the sinc function. n=4 gives similar performance to
%       a spline while n=8 is better.
%   ************ default n=8 *********
% m ... Gaussian window will be m styandard deviations down at the sinc
%       truncation point
%   ************ default m=2 *********
%
% sint ... output array of interpolated samples. There will be one column
% per input trace and the number of rows will match the length of tint.
%
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its 'AUTHOR' (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) This SOFTWARE may be used by any individual or corporation for any purpose
%    with the exception of re-selling or re-distributing the SOFTWARE.
%
% 2) The AUTHOR and CREWES must be acknowledged in any resulting publications or
%    presentations
%
% 3) This SOFTWARE is provided "as is" with no warranty of any kind
%    either expressed or implied. CREWES makes no warranties or representation
%    as to its accuracy, completeness, or fitness for any purpose. CREWES
%    is under no obligation to provide support of any kind for this SOFTWARE.
%
% 4) CREWES periodically adds, changes, improves or updates this SOFTWARE without
%    notice. New versions will be made available at www.crewes.org .
%
% 5) Use this SOFTWARE at your own risk.
%
% END TERMS OF USE LICENSE

%

%check for nan's in the interpolation times
ind=isnan(tint);
tint_original=0;
if(any(ind))
    inotnan=~isnan(tint);
    tint_original=tint;
    tint=tint(inotnan);
end

if(nargin<4)
    n=8;
end
if(nargin<5)
    m=2;
end
if(sum(abs(diff(diff(t))))>1000000*eps)
    error('input times must fall on a regular grid')
end

if(t(1)~=0)
    t1=t(1);
    t=t-t1;
    tint=tint-t1;
end

tmin=min(t);
tmax=max(t);

dt=t(2)-t(1);
nt=length(t);
[nr,nc]=size(s);
if(nr==nt)
    %the preferred case
    s2=[zeros(n,nc); s; zeros(n,nc)];%pad with zeros
    ntr=nc;
elseif(nc==nt)
    s2=[zeros(n,nr); s'; zeros(n,nr)];%pad with zeros
    ntr=nr;
else
    error('input variables t and s have incompatible sizes');
end

t2=(0:size(s2,1)-1)*dt-n*dt;
inc=50;
dt2=dt/inc;
tsinc=(-n*dt:dt2:n*dt)';
one_over_sigma=m/tsinc(end);%inverse of standard deviation of gaussian taper   
sink=sinc(tsinc/dt).*exp(-(one_over_sigma*tsinc).^2);%the sinc function table
nint=length(tint);
sint=zeros(nint,ntr);
small=100*eps;
for k=1:nint
    if(tint(k)<tmin || tint(k)>tmax)
        sint(k,:)=zeros(1,ntr);
    else
        kint=(tint(k)-t2(1))/dt+1;%fractional sample number (in s2) of the interpolation site
        if(abs(kint-round(kint))<small)
            %means the interpolation site is on an input sample
            sint(k,:)=s2(round(kint),:);
        else
            kl=floor(kint)+(-n+1:1:0);%n points to the left of tint(k)
            klsinc=(kl-1)*inc+1-round(tint(k)/dt2);%corresponding point in sink
            kr=ceil(kint)+(0:1:n-1);%n points to the right of tint(k)
            krsinc=(kr-1)*inc+1-round(tint(k)/dt2);%corresponding points in sink
            op=sink([klsinc krsinc]);%table lookup
            sint(k,:)=sum(op(:,ones(1,ntr)).*s2([kl kr],:),1);%actual interpolation
        end
    end
end
if(length(tint_original)>1)
    sint2=zeros(size(tint_original));
    sint2(inotnan)=sint;
    sint=sint2;
end

end

function y=sinc(x)
y=ones(size(x));
ind=x~=0;
xx=pi*x;
y(ind)=sin(xx(ind))./xx(ind);
end
function ind=surround(x,xtest)
% SURROUND: analyze how a vector surrounds some test points
%
% ind=surround(x,xtest)
%
% SURROUND returns a vector of indicies indicating how the vector
% x surrounds xtest which must be a scalar. If 
% isempty(ind) then xtest lies outside the range of x. Otherwise,
% ind will be the index of a point in x just greater (or less) than 
% xtest. Thus the following will be true:
%	x(ind) <= xtest < x(ind+1)
%		or
%	x(ind) >= xtest > x(ind+1)
%
% So, for if xtest is an interior point for the vector x,
% ind and ind+1 select those points in x which surround (or bracket)
% xtest. Note that x need not be monotonic. If the xtest is surrounded
% more than once by x, then ind will be a vector.
% example: x=1:10;
% >>surround(x,-1)
%    returns []
% >>surround(x,3)
%   returns 3
% now let x=[1:10 9:-1:1]
% >>surround(x,3)
%    returns 3 17
% >>surround(x,pi)
%   returns 3 16
%
% G.F. Margrave
% Jan 1994, revised Jan 95
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its 'AUTHOR' (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) This SOFTWARE may be used by any individual or corporation for any purpose
%    with the exception of re-selling or re-distributing the SOFTWARE.
%
% 2) The AUTHOR and CREWES must be acknowledged in any resulting publications or
%    presentations
%
% 3) This SOFTWARE is provided "as is" with no warranty of any kind
%    either expressed or implied. CREWES makes no warranties or representation
%    as to its accuracy, completeness, or fitness for any purpose. CREWES
%    is under no obligation to provide support of any kind for this SOFTWARE.
%
% 4) CREWES periodically adds, changes, improves or updates this SOFTWARE without
%    notice. New versions will be made available at www.crewes.org .
%
% 5) Use this SOFTWARE at your own risk.
%
% END TERMS OF USE LICENSE

	
	n=length(x);
	x1=x(1:n-1);
	x2=x(2:n);
	
	ind=find( (x1<=xtest & x2>xtest ) | ...
			(x1>=xtest & x2 < xtest ) );
end
function indicies = between(x1,x2,testpts,flag)
% BETWEEN: logical test, finds samples in vector between given bounds
%
% indicies = between(x1,x2,testpts,flag)
% indicies = between(x1,x2,testpts)
%
% returns the indicies of those points in the array testpts which lie between
% the points x1 and x2. If no testpts are found between x1 and x2 then a
% single scalar 0 (false) is returned. Flag determines the exact nature of
% the inclusion of the endpoints:
%    if flag == 0, then not endpoints are included
%    if flag == 1, then x1 is included. i.e. if a test point is precisely
%                  equal to x1,	it will be considered "between"
%    if flag == 2, then both x1 and x2 are included
%  ******** flag defaults to 0 ****
%
% Function works regardless of whether x1 < x2 or x2 < x1.
%
% by G.F. Margrave, May 1991
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its 'AUTHOR' (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) This SOFTWARE may be used by any individual or corporation for any purpose
%    with the exception of re-selling or re-distributing the SOFTWARE.
%
% 2) The AUTHOR and CREWES must be acknowledged in any resulting publications or
%    presentations
%
% 3) This SOFTWARE is provided "as is" with no warranty of any kind
%    either expressed or implied. CREWES makes no warranties or representation
%    as to its accuracy, completeness, or fitness for any purpose. CREWES
%    is under no obligation to provide support of any kind for this SOFTWARE.
%
% 4) CREWES periodically adds, changes, improves or updates this SOFTWARE without
%    notice. New versions will be made available at www.crewes.org .
%
% 5) Use this SOFTWARE at your own risk.
%
% END TERMS OF USE LICENSE

if nargin < 4, flag =0; end

if(length(x1)~=1 || length(x2)~=1)
    error('x1 and x2 must be scalars')
end

if flag == 0
    if( x1 < x2)
        indicies = find( (testpts > x1)&(testpts < x2) );
        if(isempty(indicies))
            indicies=0;
        end
        return;
    else
        indicies = find( (testpts > x2)&(testpts < x1) );
        if(isempty(indicies))
            indicies=0;
        end
        return;
    end
end
if flag == 1
    if( x1 < x2)
        indicies = find( (testpts >= x1)&(testpts < x2) );
        if(isempty(indicies))
            indicies=0;
        end
        return;
    else
        indicies = find( (testpts > x2)&(testpts <= x1) );
        if(isempty(indicies))
            indicies=0;
        end
        return;
    end
end
if flag == 2
    if( x1 < x2)
        indicies = find( (testpts >= x1)&(testpts <= x2) );
        if(isempty(indicies))
            indicies=0;
        end
        return;
    else
        indicies = find( (testpts >= x2)&(testpts <= x1) );
        if(isempty(indicies))
            indicies=0;
        end
        return;
    end
end
end