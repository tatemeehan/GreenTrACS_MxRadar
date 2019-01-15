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
    
    % Annual Accumulation Correction from Core 15 Chemistry [mwe]
    % 50 year mean is ~0.299 [mwe]
    accumulationAvg = 0.299;
    iceCoreLL = [73.592732,-47.197152];
    iceCoreXY = ll2psn(iceCoreLL(1),iceCoreLL(2));
    % Extract Radar Measurments Nearest to Ice Core.
    if isLoadGPS
        k = 100; % Number of Radar Estimates to Estimate Core Site Average
    if isCat
        iceCoreFileIx = 1;
    [~,iceCoreIx] = mink(sqrt((iceCoreXY(1)-trhd{1}(10,1:nChan:end)).^2 ...
    + (iceCoreXY(1)-trhd{1}(11,1:nChan:end)).^2),k);

    else 
        for ii = 1:nFiles
            [iceCoreDistance(ii,1:k),iceCoreIx(ii,1:k)] = mink(sqrt((iceCoreXY(1)-trhd{ii}(10,1:nChan:end)).^2 ...
    + (iceCoreXY(1)-trhd{ii}(11,1:nChan:end)).^2),k);
        end
        [~,iceCoreFileIx] = min(mean(iceCoreDistance,2));
        iceCoreIx = iceCoreIx(iceCoreFileIx,:);
    end
    else
        % If Radar is not Paired with GPS assume Ice Core is located 
        % adjacent to the first k traces of the first radar file
        iceCoreFileIx = 1;
        iceCoreIx = 1:k;        
    end
    
    % Estimate of Mean Annual Temperature    
    isLoadT2 = 1;
%     annualT = -18; % [C]     
%     annualT = -19; % [C]    
%     annualT = -20; % [C]
    annualT = -21.5; % [C]    


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
                t2Dir = 'E:\ArcticDataCenter\Data\Birkle';
                t2file = 'merra_t2_1979-2012_monthly.nc';
                [annualT] = reanalysisT2(dataDir,t2file,XY{ii});
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
            accuHL = accumulationAvg + (SnowWaterEqv{1,iceCoreFileIx}(jj) - ...
                mean(SnowWaterEqv{1,iceCoreFileIx}(iceCoreIx)));
            AverageAccumulation{ii}(jj) = accuHL;
            % Impose Herron-Langway (1980) Density Model
            [HerronLangwayDensity{ii}(:,jj),HerronLangwayAge{ii}(:,jj)] = ...
                herronLangway(StackZ,annualT,surfHL,accuHL);
            
            %%% Apply Surface Correction to Herron-Langway Model Here %%%
            radHLdepth0 = StackZ(2);
            radHLrho0 = SurfaceDensity{ii}(jj);
            [~, radHLdepthIx1] = min(abs(StackZ - dhDepth{2,ii}(jj)));
            radHLdepth1 = StackZ(radHLdepthIx1);
            radHLrho1 = dhDensity{2,ii}(jj);
            [~, radHLdepthIx2] = min(abs(StackZ - ForcingDepth{ii}(jj)));
            radHLdepth2 = StackZ(radHLdepthIx2);
            radHLrho2 = ForcingDensity{ii}(jj);
            [~, radHLdepthIx3] = min(abs(StackZ - Depth{1,ii}(jj)));
            radHLdepth3 = StackZ(radHLdepthIx3);
            radHLrho3 = Density{1,ii}(jj);
            
            % L2 Norm Regression
            G = [1,radHLdepth0; 1,radHLdepth1;1,radHLdepth2; 1,radHLdepth3];
            d = [radHLrho0;radHLrho1;radHLrho2;radHLrho3];
            m = G\d;
            surfRhoHL = m(1);
            rhoRateHL = m(2);
            
            %%% Average Radar Derived Surface Density %%%
            avgSurfacePiece = rhoRateHL.*StackZ(1:radHLdepthIx3) + surfRhoHL;
           
            %%% Causal Integration Inversion for Surface Density %%%
            G = tril(ones(length(1:radHLdepthIx3)));
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
    end