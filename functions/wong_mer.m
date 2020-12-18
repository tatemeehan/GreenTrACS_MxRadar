function [jmax, er33, SNRdb] = wong_mer(gram, winlen, sampletime, choice)
%
%function [jmax,er33, SNRdb] = wong_mer(gram, winlen, sampletim, choice);
%  Updated NOv. 24, 2010
%
% INPUT
%   gram = input seismogram
%   winlen = input energy window length in seconds
%   sampletime = in seconds (dt)
%   choice  = integer for selecting power of modified energy ratio:
%       if choice =0; er33(j) = er(j) = ordinary energy ratio;
%       if choice =1; er33(j) = er(j)^3*abs(gram(i);
%       if choice =2; er33(j)=  er(j)*abs(gram(j)) ^3 ;
%       if choice =3; er33(j) = ( er(j)*abs(gram(j)) )^3;
%       sampletime = in seconds
% OUTPUT
%   jmax = output index of time pick
%   er33 = vector showing MER ratio
%   SNRdb = SNR of enery ratio
%
% Modified by Dylan Mikesell (2/23/2012)

n = length(gram) ;

tst_len = fix(winlen/sampletime);
epsilon = 0.0001 ;
er = 0*gram + epsilon;
er33 = 0*gram + epsilon;

ap =  (gram(1)^2 +gram(2)^2)/2 ;
af =  (gram(n-1)^2 +gram(n)^2)/2 ;

for j = 1 : n
    am1= epsilon ;
    for k = j-tst_len : j
        if k <= 0
            am1 = am1 + ap ;
        else
            am1 = am1 + gram(k)*gram(k) ;
        end
    end
    
    am2= epsilon*0 ;
    for k = j : j+tst_len
        if k>=n
            am2 = am2+af ;
        else
            am2 = am2 + gram(k)*gram(k) ;
        end ;
    end
    
    am = abs(gram(j));
    ar = am2/am1 ; %  energy ratio
    er(j) = ar ;
    % modified energy ratio choices
    if choice ==0; er33(j) = ar; end
    if choice ==1; er33(j) = ar*ar*ar*am; end
    if choice ==2; er33(j)= am*am*am*ar ;  end
    if choice >=3; er33(j) = ar*ar*ar*am*am*am; end
end

amax = 0 ; jmax = 1;
for j = tst_len : n
    if amax<er33(j)
        amax = er33(j); 
        jmax= j; 
    end
end
SNRdb = 10*log10(max(er)) ;

return
    