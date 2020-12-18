function [dataOut] = SVDSfilter(data,threshold,type,BPix)
% SVDSfilter performs an economical decomposition of the GPR data. Using
% the bandpass method (Cagnoli and Ulrych 2001) coherent noise and clutter
% in the form of eigenimages are removed from the data. The first 80
% singular values are requested for a scalable SVD filter. The pass-band 
% may be specified by the user, or will be determined by PCA. Alternatively
% type may be specified for 'LP' ~ low pass or 'BP' ~ band pass.
% The variable threshold is between 0 and 1 and mutes singular values upto 
% percent PCA.
%
% Tate Meehan 12/4/2018
% Surpress Try Catch Warning
id = 'MATLAB:mir_warning_changing_try_catch';
warning('off',id);
try threshold;
catch
    threshold = 0.5;
end
try type;
catch
    % Set the Default Filter to BandPass
    type = 'BP';
end
try BPix;
catch
    % Set the Default Decomposition to 100 singular values
    BPix = [1,100];
end
% [U,S,V] = svd(data);
[U,S,V] = svds(data,BPix(2));
% Extract the Singular Values for PCA
Sigma = diag(S);           
% Principle Components of SVD
pca = Sigma./sum(Sigma);
% The cumulative PCA  value determines the filter Low Cut
% Coherent Noise is beleived to contribute between 25% - 50% of EigenImage
cumPCA = tril(ones(length(pca)))*pca;
BPix(1) = find(cumPCA >=threshold,1);
% Weight for Low Cut Eigenvalue > threshold
w = (cumPCA(BPix(1)) - threshold)./threshold;
% Kill the Largest Eigenvalues
if strcmp(type,'BP')
    if BPix(1) > 1
        S(1:BPix(1)-1,1:BPix(1)-1) = 0;
    end
    S(BPix(1),BPix(1)) = w;
elseif strcmp(type,'LP')
    % Preserve the Largest Eigenvalues
    S(BPix(1)+1:end,BPix(1)+1:end) = 0;
    S(BPix(1),BPix(1)) = w;
end
% Reconstruct the Data
dataOut = U * S * V'; 
end

