function [ cleanData, NaNno, datum ] = cleanNaN( data )
% cleanNAN removes non exitant values from a data matrix.
% NaN may arise thru data processing or lye within the raw field record. 
% To improve signal coherency across stratum, the NAN is replaced by 
% the mean of its row. 
% In signal processing this filter should be applied prior to data stacking. 
%
%   Inputs:     data - The pre-stack data marix
%
%   Output:     cleanData - The NaN replaced data
%
%   Written by Tate Meehan, Boise State University, GreenTrACS 2016

[nr,nc] = size(data);   % nr - numRow ; nc - numCol
cleanData = data;       % Allocate Data Matrix for Cleanse
NaNno = 0;              % Count Number of NaN Values Cleansed
datum = nr.*nc;
    % Clean NAN
    for mm = 1:nc
        for nn = 1:nr
            if isnan(data(nn,mm))
                cleanData(nn,mm) = mean(data(nn));
                NaNno =  NaNno + 1;
            end
        end
    end
end

