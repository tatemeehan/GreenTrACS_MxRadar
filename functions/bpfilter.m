function out = bpfilter(a,dt,fmin,fmax,btord,verbose)
%
% USAGE:    a = bpfilter(a,dt,fmin,fmax,btord,verbose)
%
% INPUT:
%   a(npts,ntrc) = data matrix
%   dt           = sample interval (s)
%   fmin,fmax    = low and high cut frequencies (Hz) (default=[1e-6,nyquist])
%   btord        = butterworth filter order (default=3)
%   verbose      = 1 if user wants information about filtering (default=0)
%
% OUTPUT:
%   a(npts,ntrc)=filtered data matrix
%
% NOTES: filtfilt function does both forward and reverse filters to ensure
% zero-phase filter.
%
% Created by Dylan Mikesell: 4 March 2013

nyqf=1/dt/2; % nyquist

% set some defaults in case user does not define anything
if nargin<2
    error('dt not set. Fix!');
end
if nargin<3
    fmin=1e-6; % set default
end
if nargin<4
    fmax=nyqf-1e-6; % set fmax to nyquist
end
if nargin<5
    btord=3; % set default butterworth order
end
if nargin<6
    verbose=0; % set default verbosity
end
if fmax/nyqf >1
    fmax = nyqf-1;
end

if verbose
    fprintf('Nyquist      = %f.\n',nyqf);
    fprintf('[fmin,fmax]  = [%f, %f].\n',fmin,fmax);
    fprintf('Filter order = %d.\n',btord);
end

% build filter
[z,q]=butter(btord,[fmin/nyqf fmax/nyqf],'bandpass');
[nsamp,ntrc]=size(a);
out = zeros(nsamp,ntrc);
% if isempty(gcp('nocreate')) % serial
%     for ii=1:ntrc
% %         tmp=demean(a(:,ii));
%         tmp=a(:,ii)-mean(a(:,ii));
%         tmp=detrend(tmp);
%         w = tukeywin(numel(tmp),0.05); % taper
%         a(:,ii)=filtfilt(z,q,tmp.*w);
%     end
% else % or parallel
%     parfor ii=1:ntrc
for ii = 1:ntrc
    tmp = a(:,ii);
    tmp=a(:,ii)-mean(a(:,ii));
    w = tukeywin(numel(tmp),0.05); % taper
    tmp=filtfilt(z,q,tmp.*w);
    out(:,ii) = tmp(1:nsamp);
end
% end

return