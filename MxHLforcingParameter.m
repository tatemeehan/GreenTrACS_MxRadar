%% MxHL Forcing Parameterization
% This Routine Estimates the Surface Density and Forcing Density for the 
% MxHL model and the uncertainties of these paramters. 

if ~isSWEDISH || isLoadHVA
    % Determine Maximum File Size
    nTrace = zeros(1,nFiles);
    for ii = 1:nFiles
        nTrace(ii) = size(Radar{1,ii},2);
    end
    % Maxiumum Traces to Loop Over
    nTrace = max(nTrace);
    
    % Smoothing Window Length
    smoothR = 251;
end

if ~isLoadHVA || ~isLoadMxHL
% Allocation
looper = 1:nTrace;
SurfaceDensity = cell(nFiles,1);
SurfaceDensityVar = cell(nFiles,1);
SurfaceDensification = cell(nFiles,1);
SurfaceDensificationVar = cell(nFiles,1);
ForcingDensity = cell(nFiles,1);
ForcingDensityVar = cell(nFiles,1);
ForcingDepth = cell(nFiles,1);
ForcingDepthVar = cell(nFiles,1);

for ii = 1:nFiles
    MCsamples = 250;
    % Allocation
    fDensity = zeros(length(looper),1);
    fDensityVar = zeros(length(looper),1);
    fDepth = zeros(length(looper),1);
    fDepthVar = zeros(length(looper),1);
    surfVar = zeros(length(looper),1);
    rateVar = zeros(length(looper),1);
    surfDensity = zeros(length(looper),1);
    densityRate = zeros(length(looper),1);
    
    % MC BootStrapping Regression for Forcing Density and Surface Density
    parfor (jj = looper, nWorkers)

        % Ensure Real Indicies
        realIx = find(real([xToRef{1,jj,ii}]));
        % Create Indicies for Sampling
        dirIx = 1:length([xRhoDir{2,jj,ii}]);
        refIx = 1:length([xRhoRef{1,jj,ii}(realIx)]);
        
        % Monte Carlo Sampling Depth and Density
        dirSample = randsample(dirIx,MCsamples,'true');
        refSample = randsample(refIx,MCsamples,'true');
        
        % dh = 2 Surface Wave Horizon
        xStrapDepthDir = xDepthDir{2,jj,ii}(dirSample);
        xStrapRhoDir = xRhoDir{2,jj,ii}(dirSample);

        % rh = 1 Primary Reflection Horizon
        xStrapDepthRef = xDepth{1,jj,ii}(refSample);
        xStrapRhoRef = xRhoRef{1,jj,ii}(refSample);
        
        % The Forcing Density of Herron-Langway (1980) is to Good Approximation the
        % the mean of the LMO and NMO derived density
        fDensity(jj) = mean(mean([xStrapRhoDir,xStrapRhoRef],2));
        fDensityVar(jj) = var(mean([xStrapRhoDir,xStrapRhoRef],2));
        fDepth(jj) =  mean(mean([xStrapDepthDir,xStrapDepthRef],2));
        fDepthVar(jj) = var(mean([xStrapDepthDir,xStrapDepthRef],2));
        
        % Allocation
        surfRho = zeros(MCsamples,1);
        rhoRate = zeros(MCsamples,1);
        for kk = 1:MCsamples
            % L2 Norm Regression
            G = [1,xStrapDepthDir(kk);1,xStrapDepthRef(kk)];
            d = [xStrapRhoDir(kk);xStrapRhoRef(kk)];
            m = G\d;
            surfRho(kk) = m(1);
            rhoRate(kk) = m(2);
        end
        % Estimate BootStrapped Sample Variance
        surfVar(jj) = var(surfRho);
        rateVar(jj) = var(rhoRate);
        surfDensity(jj) = mean(surfRho);
        densityRate(jj) = mean(rhoRate);
    end
    SurfaceDensity{ii} = nonParametricSmooth( 1:length(surfDensity),...
        surfDensity,1:length(surfDensity),smoothR);
    SurfaceDensityVar{ii} = nonParametricSmooth( 1:length(surfVar),...
        surfVar,1:length(surfVar),smoothR);
    SurfaceDensification{ii} = nonParametricSmooth( 1:length(densityRate),...
        densityRate,1:length(densityRate),smoothR);
    SurfaceDensificationVar{ii} = nonParametricSmooth( 1:length(rateVar),...
        rateVar,1:length(rateVar),smoothR);

% The Forcing Density of Herron-Langway (1980) is to Good Approximation the
% the mean of the LMO and NMO derived density
ForcingDensity{ii} = nonParametricSmooth( 1:length(fDensity),...
        fDensity,1:length(fDensity),smoothR);
ForcingDensityVar{ii} = nonParametricSmooth( 1:length(fDensityVar),...
        fDensityVar,1:length(fDensityVar),smoothR);
ForcingDepth{ii} = nonParametricSmooth( 1:length(fDepth),...
        fDepth,1:length(fDepth),smoothR);
ForcingDepthVar{ii} = nonParametricSmooth( 1:length(fDepthVar),...
        fDepthVar,1:length(fDepthVar),smoothR);
    
end
end