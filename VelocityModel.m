%% Velocity Model produces the MxHL firn desnity model.
% The Herron and Langway (1980) model requires input properties for average
% annual snow accumulation, average annual temperature, and the average 
% 1 - 2 m snow density. GreenTrACS ice cores are chemically dated by a 
% continuous melting process (Osterberg et al., 2006; Graeter et al. 2017),
% and provide the known average accumulation through the time period of
% modern reanalyis models. The radar estimates are tranfered to this
% calibration and further estimate deviations from the (~1980 - 2017 mean).
% Average annual temperature is estimated from GreenTrACS reanalysis data 
% available via the Arctic Data Center: 
% Sean Birkel. 2018. Greenland surface mass balance derived from climate 
% reanalysis models, 1979-2017. Arctic Data Center. doi:10.18739/A2D21RH75.
%
% Tate Meehan, Boise State University, NASA ISGC 2019
% This Routine Build the 2-D EM Velocity Model 
    % Allocation
    AverageAccumulation = cell(nFiles,1);
    StackingVelocity = cell(nFiles,1);
    StackingDensity  = cell(nFiles,1);
    StackingTime = cell(nFiles,1);
    HerronLangwayDensity = cell(nFiles,1);
    HerronLangwayAge = cell(nFiles,1);
    Traverse = cell(nFiles,1);
    TraverseX = cell(nFiles,1);
    StackDepth = cell(nFiles,1);
    XY = cell(nFiles,1);
    
    % Estimate of Mean Annual Temperature    
    isLoadT2 = 1;
    if ~isLoadT2
      annualT = -25.5; %[C]
    end

    % Array of Model Depths
    maxDepth = 30.001;
    nDepth = 301;
    StackZ = (linspace(0.001,maxDepth,nDepth))';
    % Interval Thickness
    StackH = diff(StackZ);StackH = [StackZ(1);StackH];
    
    for ii = 1:nFiles
        if isLoadGPS
            % Coordinate Positions
            XY{ii} = [trhd{ii}(10,1:nChan:end);trhd{ii}(11,1:nChan:end)];
            Traverse{ii} = trhd{ii}(2,1:nChan:end);
            % Distance Grid
            TraverseX{ii} = ones(length(StackZ),1)*Traverse{ii};
            % Depth Grid
            StackDepth{ii} = StackZ*ones(1,length(Traverse{ii}));
            
            if isLoadT2
                t2Dir = '/home/tatemeehan/GreenTracs2017/ClimateData/';
%                 t2Dir = 'D:\ArcticDataCenter\Data\Birkle';
                t2file = 'merra_t2_1979-2012_monthly.nc';
                [annualT] = reanalysisT2(t2Dir,t2file,XY{ii});
            end
        else
            % Array of Approximate Distance
            Traverse{ii} = linspace(0,TraverseDistance(ii),length(dhDensity{2,1}));
            % Distance Grid
            TraverseX{ii} = ones(length(StackZ),1)*Traverse{ii};
            % Depth Grid
            StackDepth{ii} = StackZ*ones(1,length(Traverse{ii}));
        end
        
        for jj = 1:length(ForcingDensity{ii})
            % Extract Radar Forcing for Herron-Langway Model
            surfHL = ForcingDensity{ii}(jj);
            % radar estimates are tranfered to ice core accumulation 
            accuHL = coreAccumulation(1) + (SnowWaterEqv{1,iceCoreFileIx}(jj) - ...
                mean(SnowWaterEqv{1,iceCoreFileIx}(iceCoreIx)));
            AverageAccumulation{ii}(jj) = accuHL;
            % Impose Herron-Langway (1980) Density Model
            [HerronLangwayDensity{ii}(:,jj),HerronLangwayAge{ii}(:,jj)] = ...
                herronLangway(StackZ,annualT(jj),surfHL,accuHL);
            % Compute Cumulative Average
            avgHerronLangwayDensity = movmean(HerronLangwayDensity{ii}(:,jj),...
                [length(HerronLangwayDensity{ii}(:,jj))-1 0]);
            %%% Apply Surface Correction to Herron-Langway Model Here %%%
            radHLdepth0 = StackZ(2);
            radHLrho0 = SurfaceDensity{ii}(jj);
            [~, radHLdepthIx1] = min(abs(StackZ - ForcingDepth{ii}(jj)));
            radHLdepth1 = StackZ(radHLdepthIx1);
            radHLrho1 = avgHerronLangwayDensity(radHLdepthIx1);
            
            % L2 Regression
            G = [1,radHLdepth0; 1,radHLdepth1];
            d = [radHLrho0;radHLrho1];
            m = G\d;
            surfRhoHL = m(1);
            rhoRateHL = m(2);
            
            %%% Average Radar Derived Surface Density %%%
            avgSurfacePiece = rhoRateHL.*StackZ(1:radHLdepthIx1) + surfRhoHL;
           
            %%% Causal Integration Inversion for Surface Density %%%
            G = tril(ones(length(1:radHLdepthIx1)));
            d = avgSurfacePiece.*sum(G,2);
            m = G\d;
            SurfacePiece = m;
            
            % Find Intersection
            [~, PieceIx] = min(abs(HerronLangwayDensity{ii}(1:length(SurfacePiece),jj)-SurfacePiece));

            % Append Surface Piece
            HerronLangwayDensity{ii}(1:PieceIx-1,jj) = SurfacePiece(1:PieceIx-1);
            
            % Compute Weights for Cumulative Average Density
            stackW = StackZ.^(-1).*tril(ones(length(StackH),1)*StackH');

            % Cumulative Average Density as a Proxy for Stacking Velocity
            StackingDensity{ii}(:,jj) = stackW*HerronLangwayDensity{ii}(:,jj);
            
            % Estimate Stacking Velocity in Time
            StackingVelocity{ii}(:,jj) = DryCrimVRMS(StackingDensity{ii}(:,jj));
            StackingTime{ii}(:,jj) = 2.*StackZ./StackingVelocity{ii}(:,jj);
                      
        end
        % Apply Correction to Interannual Isochrone Depths
        if isGreenTracsFirnCore
            % Currently works for 1 Firn Core in Analysis
            % Shift from Calendar Years to Ages
            tmp1 = abs(GTCdepthAge(1,2)-GTCdepthAge(:,2)); 
            % Extract Age-Depth Model Near Firn Core and Average
            tmp2 = mean(HerronLangwayAge{ii}(:,iceCoreIx(ii,:)),2);
            % Compute Residual (Observed - Estimated)
            isochroneResidual(:,ii) = tmp1-tmp2;
            % Update Model Positive is Downwars so Residual should be Added
            HerronLangwayAge{ii} = HerronLangwayAge{ii} + isochroneResidual(:,ii);
            % Remove any Numerical Bias
            HerronLangwayAge{ii}(1,:) = zeros(1,size(HerronLangwayAge{ii}(1,:),2));
        end
    end