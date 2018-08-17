%% Velocity Model
% This Routine Build the 2-D Velocity Model
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
    
    % Annual Accumulation Correction from Pit 15W [.Percent]
    % This is the depth Bias Correction Factor
%     normSWE = 0.0990;
    % Annual Accumulation Correction from Core 15 Chemistry [mwe]
    % 50 year mean is ~0.299 [mwe]
    accumulationAvg = 0.299;


    % Estimate of Mean Annual Temperature
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
        % Array of Approximate Distance
        Traverse{ii} = linspace(0,TraverseDistance(ii),length(dhDensity{2,1}));
        % Distance Grid
        TraverseX{ii} = ones(length(StackZ),1)*Traverse{ii};
        % Depth Grid
        StackDepth{ii} = StackZ*ones(1,length(Traverse{ii}));
        for jj = 1:length(ForcingDensity{ii})
            % Extract Radar Forcing for Herron-Langway Model
            surfHL = ForcingDensity{ii}(jj);
%             accuHL = SnowWaterEqv{1,ii}(jj);% - (normSWE.*SnowWaterEqv{1,ii}(jj));
            
            accuHL = accumulationAvg + (SnowWaterEqv{1,ii}(jj) - mean(SnowWaterEqv{1,ii}(1:100)));
            AverageAccumulation{ii}(jj) = accuHL;
            % Impose Herron-Langway (1980) Density Model
            [HerronLangwayDensity{ii}(:,jj),HerronLangwayAge{ii}(:,jj)] = herronLangway(StackZ,annualT,surfHL,accuHL);
            
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

            % Cumulative Average Density is a Proxy for Stacking Velocity
            StackingDensity{ii}(:,jj) = stackW*HerronLangwayDensity{ii}(:,jj);
            
            % Estimate Stacking Velocity in Time
            StackingVelocity{ii}(:,jj) = DryCrimVRMS(StackingDensity{ii}(:,jj));
            StackingTime{ii}(:,jj) = 2.*StackZ./StackingVelocity{ii}(:,jj);
                      
        end
    end