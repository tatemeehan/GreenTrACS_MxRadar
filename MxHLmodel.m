%% MxHL Model
% Builds the MxHL Model, Age-Depth, and Accumulation Models to match the
% size of the Radar Image.

    % Interpolate Density Model to Resolution of Radar Depth Section
    AgeModel = cell(nFiles,1);
    Ages = cell(nFiles,1);
    DensityModel = cell(nFiles,1);
    DensityAnomalyModel = cell(nFiles,1);
    AvgDensityModel = cell(nFiles,1);
    MeanDensityDeviation = cell(nFiles,1);
    SurfaceDensityDeviation = cell(nFiles,1);
    DepthAge = cell(nFiles,1);
    AvgAccumulation = cell(nFiles,1);
    MxAccumulationRate = cell(nFiles,1);
    MxHLAccumulationRate = cell(nFiles,1);
    ix550 = zeros(1,size(DensityModel{ii},2));
    depth550 = zeros(1,size(DensityModel{ii},2));
    for ii = 1:nFiles
        % Resize Depth-Density Model
        DensityModel{ii} = griddata(TraverseX{ii},StackDepth{ii},...
            HerronLangwayDensity{ii},xStack{ii},DepthMatrix{ii},'linear');
        % Fill in NaN with Nearest Extrapolation
        DensityModel{ii}(1,:) = DensityModel{ii}(2,:);
        
        % Resize Depth-AvgDensity Model
        AvgDensityModel{ii} = griddata(TraverseX{ii},StackDepth{ii},...
            StackingDensity{ii},xStack{ii},DepthMatrix{ii},'linear');
        % Fill in NaN with Nearest Extrapolation
        AvgDensityModel{ii}(1,:) = AvgDensityModel{ii}(2,:);        
        
        % Resize Depth-Age Model
        AgeModel{ii} = griddata(TraverseX{ii},StackDepth{ii},...
            HerronLangwayAge{ii},xStack{ii},DepthMatrix{ii},'linear');
        % Fill in NaN with Zero Age
        AgeModel{ii}(1,:) = 0;
        % Determine Critical Depth 
        for kk = 1:size(DensityModel{ii},2)
            ix550(kk) = find(DensityModel{ii}(:,kk)>.550,1);
            depth550(kk) = DepthMatrix{ii}(ix550(kk),kk);
        end
        
        % Density Anomaly from the Mean
        DensityAnomalyModel{ii} = DensityModel{ii} - mean(DensityModel{ii},2);

        MeanDensityDeviation{ii} = mean(DensityAnomalyModel{ii},1);
        SurfaceDensityDeviation{ii} = DensityAnomalyModel{ii}(1,:);

        % Calculate Accumulation Rate from Maximum Isochrone to Surface
        
        MaxAge = floor(min(AgeModel{ii}(end,:)));
        Ages{ii} = [1,3:3:MaxAge];
        Isochrones = abs(Year(ii) - Ages{ii});
        for kk = 1:length(Isochrones)
        [~, AgeIx] = min(abs(AgeModel{ii}-Ages{ii}(kk)));
        AgeIx = sub2ind(size(DepthMatrix{ii}),AgeIx,1:length(AgeIx));
        DepthAge{ii}(:,kk) = DepthMatrix{ii}(AgeIx');
        AvgAccumulation{ii}(:,kk) = AvgDensityModel{ii}(AgeIx').*DepthAge{ii}(:,kk);
        end
        % This Calculation is Redundant Due Circularity of Methodology
        MxHLAccumulationRate{ii} = AvgAccumulation{ii}(:,find(Ages{ii} == MaxAge))./MaxAge;
        % MxRadar Evaluated Contemporary Accumulation Rate
        % This result is Bias Corrected (Depth and Density) is Final
        MxAccumulationRate{ii} = SnowWaterEqv{1,ii};% - (normSWE.*SnowWaterEqv{1,ii});
        % Uncertainties Are Developed During HVA Bootstrapping
        Isochrones = DepthAge;
    end