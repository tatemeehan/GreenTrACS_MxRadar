%% MxHL Forcing Parameterization
% This Routine Estimates the Surface Density and Forcing Density for the 
% MxHL modeland the uncertainties of these paramters. 

if ~isSWEDISH
    % Determine Maximum File Size
    nTrace = zeros(1,nFiles);
    for ii = 1:nFiles
        nTrace(ii) = size(Radar{1,ii},2);
    end
    % Maxiumum Traces to Loop Over
    nTrace = max(nTrace);
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
    surfVar = zeros(length(looper),1);
    rateVar = zeros(length(looper),1);
    surfDensity = zeros(length(looper),1);
    densityRate = zeros(length(looper),1);
    
    % MC BootStrapping Regression for Surface Density Extrapolation
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
        
        % Allocation
        surfRho = zeros(MCsamples,1);
        rhoRate = zeros(MCsamples,1);
        for kk = 1:250
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
% Average of the two Radar Surface Estimates.
ForcingDensity{ii} = mean([dhDensity{2,ii},Density{1,ii}],2);
ForcingDensityVar{ii} = mean([dhDensityVar{2,ii},DensityVar{1,ii}],2);
ForcingDepth{ii} = mean([dhDepth{2,ii},Depth{1,ii}],2);
ForcingDepthVar{ii} = mean([dhDepthVar{2,ii},DepthVar{1,ii}],2);
    
end
end