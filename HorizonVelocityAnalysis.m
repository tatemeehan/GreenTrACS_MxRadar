%% Horizon Velocity Analysis
% This Routine Extracts the Travel-Time Picks and computes the Monte-Carlo
% Bootstrapped Linear Regression for EM Velocity and Intercept Time.
% These Parameters are Converted to Depth, Density, and Accumualtion with
% their respective variances in this routine.

% Surpress Parfor Temporary Variable Warning
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');

% Extract AirWave Picks and Perform Residual Subtraction Velocity Analysis
%     parfor (ii = 1:nFiles, nWorkers - (nWorkers-nFiles)) 
    for ii = 1:nFiles
          looper = 1:nTrace;
        % Concatenate Air Wave Picks
        for dh = 1:nDirectHorizon
            % Concatenate Primary Reflection Picks for Horizon hh
            GatherDirectPicks{dh,ii} = cat(2,DirectFBpick{:,dh,ii});
            DirectPicks = GatherDirectPicks{dh,ii};
            % Air Wave Velocity Analysis
            if dh  == 1
%                 for jj = looper
                 parfor (jj = looper, nWorkers)
                    
                    % Extract Picks
                    AirPick = DirectPicks(jj,:);
                    
                    for kk = 1:250 % 250 Random Draws
                        % Cross-Validation for Surface Velocity Estimation 7-28-17
                        nCut = randsample([0,1,2],1);
                        cutChan = randsample(liveChan,nCut);
                        xvalChan = liveChan;
                        cutIx = find(ismember(liveChan,cutChan));
                        xvalChan(cutIx) = [];
                        xvalIx = find(ismember(liveChan,xvalChan));
                        
                        xvalAirPick = AirPick(xvalIx);
                        xvalOffset = offsetArray(xvalIx);
                        
                        % Compute Air Wave Arrival Velocity and Intercept Time
                        % IRLS Scheme
                        if isL1Air
                            [xVdir{dh,jj,ii}(kk,1), xToDir{dh,jj,ii}(kk,1)] = DirectWaveIrls(xvalOffset,xvalAirPick);
                        end
                        % OLS Scheme
                        if isL2Air
                            [xVdir{dh,jj,ii}(kk,1), xToDir{dh,jj,ii}(kk,1)] = DirectWave(xvalOffset,xvalAirPick);
                        end
                       
                    end
                    
                    % Residual Error Analysis
                    VoAir = mean([xVdir{dh,jj,ii}]);
                    tmpTo = mean([xToDir{dh,jj,ii}]);
                    
                    % Find Velocity Residual for Additional Trace Shifts & Calulate Residual Time
                    deltaT{ii,jj} = (AirPick-tmpTo) - (offsetArray./0.299);
                    % Apply Residual Time Shifts
                    ResAirPick = AirPick - deltaT{ii,jj};

                    % Recalculate Airwave Moveout Velocity
                    if isL1Air
                        [VoDir{dh,jj,ii}, toDir{dh,jj,ii}] = DirectWaveIrls(offsetArray,ResAirPick);
                    end
                    if isL2Air
                        [VoDir{dh,jj,ii}, toDir{dh,jj,ii}] = DirectWave(offsetArray,ResAirPick);
                    end

                    AirTo{ii,jj} = toDir{dh,jj,ii};

                    % Estimate Wavelet Depth
                    xDepthDir{dh,jj,ii} = 0;
                    % Estimate Direct Wave Density
                    xRhoDir{dh,jj,ii} = 0;
                    % Error Analysis
                    VoVarDir{dh,jj,ii} = var([xVdir{dh,jj,ii}]);
                    toVarDir{dh,jj,ii} = var([xToDir{dh,jj,ii}]);
                    DepthDir{dh,jj,ii} = VoDir{dh,jj,ii}.*toDir{dh,jj,ii};
                    RhoDir{dh,jj,ii} = DryCrim(VoDir{dh,jj,ii});
                    CovZP = zeros(2,2);
                    DepthVarDir{dh,jj,ii} = CovZP(1,1);
                    RhoDirVar{dh,jj,ii} = CovZP(2,2);
                    CovDepthRhoDir{dh,jj,ii} = CovZP(1,2);
                end
            else
                % Cross-Validation for Surface Velocity Estimation
%                 for jj = looper
                 parfor (jj = looper, nWorkers)
                    % Multiple Shot Gathers in Population
                    isManyShots = 1;
                    if isManyShots
                        shotRange = 5;
                        ranger = sqrt(([1:nTrace]-jj).^2); % Compute Distance
                        getIx = find(ranger<=shotRange); % Find Nearby Picks
                        DirPick = DirectPicks(getIx,:) - vertcat(deltaT{ii,getIx}); % Residual Static Corection
                        pickPool = size(DirPick,1);
                        xvalPool = 1:pickPool;
                        % Cross-Validation for Surface Velocity Estimation
                        for kk = 1:1250 % 1250 Random Draws
                            nCut = randsample([0,1],1);
                            cutChan = randsample(liveChan,nCut);
                            xvalChan = liveChan;
                            cutIx = find(ismember(liveChan,cutChan));
                            xvalChan(cutIx) = [];
                            xvalIx = find(ismember(liveChan,xvalChan));
                            pickPool = length(xvalIx);
                            xvalBin = randsample(xvalPool,pickPool,'true');
                            xvalDirPick = ...
                                DirPick(sub2ind(size(DirPick),xvalBin,xvalIx));
                            xvalOffset = offsetArray(xvalIx);
                            
                            % Compute Surface Wave Arrival Velocity and Intercept Time
                            if isL2LMO
                            % OLS Scheme
                            [xVdir{dh,jj,ii}(kk,1), xToDir{dh,jj,ii}(kk,1)]...
                              = DirectWave(xvalOffset,xvalDirPick);
                            end
                            % IRLS Scheme
                            if isL1LMO
                            [xVdir{dh,jj,ii}(kk,1), xToDir{dh,jj,ii}(kk,1)]...
                                = DirectWaveIrls(xvalOffset,xvalDirPick);
                            end
                        end
                    end
                    % Single Shot Gather in Population
                    isSingleShot = 0;
                    if isSingleShot
                        % Re-Extract Direct Wave Picks after Residual Shifts
                        DirPick = DirectPicks(jj,:) - deltaT{ii,jj}; % Residual Static Corection
                        for kk = 1:250 % 250 Random Draws
                            % Cross-Validation for Surface Velocity Estimation
                            nCut = randsample([0,1,2],1);
                            cutChan = randsample(liveChan,nCut);
                            xvalChan = liveChan;
                            cutIx = find(ismember(liveChan,cutChan));
                            xvalChan(cutIx) = [];
                            xvalIx = find(ismember(liveChan,xvalChan));
                            
                            xvalDirPick = DirPick(xvalIx);
                            xvalOffset = offsetArray(xvalIx);
                            
                            % Compute Air Wave Arrival Velocity and Intercept Time
                            if isL2LMO
                            % OLS Scheme
                            [xVdir{dh,jj,ii}(kk,1), xToDir{dh,jj,ii}(kk,1)]...
                              = DirectWave(xvalOffset,xvalDirPick);
                            end
                            % IRLS Scheme
                            if isL1LMO
                            [xVdir{dh,jj,ii}(kk,1), xToDir{dh,jj,ii}(kk,1)]...
                                = DirectWaveIrls(xvalOffset,xvalDirPick);
                            end
                        end
                    end
                    % Correct Residual Direct Arrival Intercept Time
                    if dh == 2
                        ResidualToDir = AirTo{ii,jj} - mean(xToDir{dh,jj,ii});
                        AirTo{ii,jj} = AirTo{ii,jj} - ResidualToDir;
                    end
                    % Apply To Static Shift Post Residual Correction
                    xToDir{dh,jj,ii} = xToDir{dh,jj,ii} - AirTo{ii,jj};
                    % Apply Velocity Bias Correction Factor 1.72%
                    xVdir{dh,jj,ii} = xVdir{dh,jj,ii}%; + velocityBias;
                    % Estimate Wavelet Depth
                    % Sounding Depth at One Wavelength
                    waveLength = (xVdir{dh,jj,ii}./f0GHz);
                        % L1 Norm Forces Variance of To Positive Downward
                    xDepthDir{dh,jj,ii} = xVdir{dh,jj,ii}.*abs(xToDir{dh,jj,ii})...
                        + waveLength;
                    % Estimate Direct Wave Density
                    xRhoDir{dh,jj,ii} = DryCrim([xVdir{dh,jj,ii}]);
                    % Error Analysis
                    VoDir{dh,jj,ii} = mean([xVdir{dh,jj,ii}]);
                    VoVarDir{dh,jj,ii} = var([xVdir{dh,jj,ii}]);
                    toDir{dh,jj,ii} = mean(abs([xToDir{dh,jj,ii}]));
                    toVarDir{dh,jj,ii} = var([xToDir{dh,jj,ii}]);
                    DepthDir{dh,jj,ii} = mean([xDepthDir{dh,jj,ii}]);
                    RhoDir{dh,jj,ii} = mean([xRhoDir{dh,jj,ii}]);
                    CovZP = cov([xDepthDir{dh,jj,ii}],[xRhoDir{dh,jj,ii}]);
                    DepthVarDir{dh,jj,ii} = CovZP(1,1);
                    RhoDirVar{dh,jj,ii} = CovZP(2,2);
                    CovDepthRhoDir{dh,jj,ii} = CovZP(1,2);
                 end
                % Apply Residual Corrections to AirWave Arrival Time
                if dh == 2
                    toDir(1,:,ii) = AirTo(ii,:);
                    DirectTo{1,ii} = [toDir{1,:,ii}] - [AirTo{ii,:}];
                end
            end
            
            % Concatenate Direct Travel Time
            DirectTo{dh,ii} = [toDir{dh,:,ii}];
            DirectToVar{dh,ii} = [toVarDir{dh,:,ii}];
            % Concatenate Direct Velocity
            DirectVelocity{dh,ii} = [VoDir{dh,:,ii}];
            DirectVelocityVar{dh,ii} = [VoVarDir{dh,:,ii}];
            % Concatenate Direct Wave Depth
            DirectDepth{dh,ii} = [DepthDir{dh,:,ii}];
            DirectDepthVar{dh,ii} = [DepthVarDir{dh,:,ii}];
            % Concatenate Reflection Density
            DirectDensity{dh,ii} = [RhoDir{dh,:,ii}];
            DirectDensityVar{dh,ii} = [RhoDirVar{dh,:,ii}];
            % Concatenate Covariance
            CovarianceDepthDensityDirect{dh,ii} = [CovDepthRhoDir{dh,:,ii}];
            
            % Estimte DirectSWE
            DirectSWE{dh,ii} = DirectDepth{dh,ii}.*DirectDensity{dh,ii};
            % Error Propagation Equation
            DirectSWEvar{dh,ii} = DirectDepthVar{dh,ii}.*DirectDensity{dh,ii}.^2 ...
                + DirectDensityVar{dh,ii}.*DirectDepth{dh,ii}.^2 ...
                + 2.*DirectDepth{dh,ii}.*DirectDensity{dh,ii}.*CovarianceDepthDensityDirect{dh,ii};
            
            % Smooth Esitmates
            smoothR = 251;
            dhTWT{dh,ii} = nonParametricSmooth( 1:length(DirectTo{dh,ii}),...
                DirectTo{dh,ii},1:length(DirectTo{dh,ii}),smoothR);
            dhTWTvar{dh,ii} = nonParametricSmooth( 1:length(DirectToVar{dh,ii}),...
                DirectToVar{dh,ii},1:length(DirectToVar{dh,ii}),smoothR);
            dhSnowWaterEqv{dh,ii} = nonParametricSmooth( 1:length(DirectSWE{dh,ii}),DirectSWE{dh,ii},...
                1:length(DirectSWE{dh,ii}),smoothR);
            dhSnowWaterEqvVar{dh,ii} = nonParametricSmooth( 1:length(DirectSWEvar{dh,ii}),...
                DirectSWEvar{dh,ii},1:length(DirectSWEvar{dh,ii}),smoothR);
            dhDensity{dh,ii} = nonParametricSmooth( 1:length(DirectDensity{dh,ii}),...
                DirectDensity{dh,ii},1:length(DirectDensity{dh,ii}),smoothR);
            dhDensityVar{dh,ii} = nonParametricSmooth( 1:length(DirectDensityVar{dh,ii}),...
                DirectDensityVar{dh,ii},1:length(DirectDensityVar{dh,ii}),smoothR);
            dhDepth{dh,ii} = nonParametricSmooth( 1:length(DirectDepth{dh,ii}),...
                DirectDepth{dh,ii},1:length(DirectDepth{dh,ii}),smoothR);
            dhDepthVar{dh,ii} = nonParametricSmooth( 1:length(DirectDepthVar{dh,ii}),...
                DirectDepthVar{dh,ii},1:length(DirectDepthVar{dh,ii}),smoothR);
        end
    end

% Primary Reflection Velocity Analysis               
    for ii = 1:nFiles
        looper = 1:nTrace;
        for rh = 1:nReflectionHorizon
        % Concatenate Primary Reflection Picks for Horizon rh
        GatherReflectionPicks{rh,ii} = cat(2,ReflectionFBpick{:,rh,ii});
        Reflection = GatherReflectionPicks{rh,ii};
%         for jj = looper
        parfor (jj = looper, nWorkers)
            % Jackknife Simulation for Reflection Velocity Estimation
            % Multiple Shot Gathers in Population
            isManyShots = 1;
            if isManyShots
                shotRange = 5;
                ranger = sqrt(([1:nTrace]-jj).^2); % Compute Distance   
                getIx = find(ranger<=shotRange); % Find Nearby Picks
                RefPick = Reflection(getIx,:) - vertcat(deltaT{ii,getIx}); % Residual Static Corection
                RefPick = RefPick - vertcat(AirTo{ii,getIx})*ones(1,nChan); % Time-Zero Static Shift
                pickPool = size(RefPick,1);
                xvalPool = 1:pickPool;
            end
            % Single Shot Gather in Population
            isSingleShot = 0;
            if isSingleShot
                RefPick = Reflection(jj,:) - deltaT{ii,jj}; % Residual
                RefPick = RefPick - AirTo{ii,jj};           % Bulk Shift
            end
                % Jackknife Simulation for Reflection Velocity Estimation
                if isManyShots
                    for kk = 1:1250 % 1250 Random Draws
                    nCut = randsample([0,1,2],1);
                    cutChan = randsample(liveChan,nCut);
                    xvalChan = liveChan;
                    cutIx = find(ismember(liveChan,cutChan));
                    xvalChan(cutIx) = [];
                    xvalIx = find(ismember(liveChan,xvalChan));
                    pickPool = length(xvalIx);
                    xvalBin = randsample(xvalPool,pickPool,'true');
                    xvalRefPick = ...
                        RefPick(sub2ind(size(RefPick),xvalBin,xvalIx));
                    xvalOffset = offsetArray(xvalIx);
                    
                    % Compute Reflected Arrival Velocity and Intercept Time
                    % IRLS Scheme
                    if isL1NMO
                        [xVref{rh,jj,ii}(kk,1), xToRef{rh,jj,ii}(kk,1), xDepth{rh,jj,ii}(kk,1)] ...
                            = VrmsIrls(xvalOffset,xvalRefPick);
                    end
                    % OLS Scheme
                    if isL2NMO
                        [xVref{rh,jj,ii}(kk,1), xToRef{rh,jj,ii}(kk,1), xDepth{rh,jj,ii}(kk,1)] ...
                            = Vrms(xvalOffset,xvalRefPick);
                    end
                    end
                end
                if isSingleShot
                    for kk = 1:250 % 250 Random Draws
                        nCut = randsample([0,1,2],1);
                        cutChan = randsample(liveChan,nCut);
                        xvalChan = liveChan;
                        cutIx = find(ismember(liveChan,cutChan));
                        xvalChan(cutIx) = [];
                        xvalIx = find(ismember(liveChan,xvalChan));
                        
                        xvalRefPick = RefPick(xvalIx);
                        xvalOffset = offsetArray(xvalIx);                    
                    % Compute Reflected Arrival Velocity and Intercept Time
                    % IRLS Scheme
                    if isL1NMO
                        [xVref{rh,jj,ii}(kk,1), xToRef{rh,jj,ii}(kk,1), xDepth{rh,jj,ii}(kk,1)] ...
                            = VrmsIrls(xvalOffset,xvalRefPick);
                    end
                    % OLS Scheme
                    if isL2NMO
                        [xVref{rh,jj,ii}(kk,1), xToRef{rh,jj,ii}(kk,1), xDepth{rh,jj,ii}(kk,1)] ...
                            = Vrms(xvalOffset,xvalRefPick);
                    end
                    end
                end
            % Apply Velocity Bias Correction Factor % 1.72
            xVref{rh,jj,ii} = xVref{rh,jj,ii};% + velocityBias;
            % Create Bootstrapped Population of Reflection Density
            xRhoRef{rh,jj,ii} = DryCrim([xVref{rh,jj,ii}]);
            
            % Error Analysis
            realIx = find(real([xToRef{rh,jj,ii}]));
            VoRef{rh,jj,ii} = mean([xVref{rh,jj,ii}(realIx)]);
            VoVarRef{rh,jj,ii} = var([xVref{rh,jj,ii}(realIx)]);
            toRef{rh,jj,ii} = mean([xToRef{rh,jj,ii}(realIx)]);
            toVarRef{rh,jj,ii} = var([xToRef{rh,jj,ii}(realIx)]);
            HorizonDepth{rh,jj,ii} = mean([xDepth{rh,jj,ii}(realIx)]);
            RhoRef{rh,jj,ii} = mean([xRhoRef{rh,jj,ii}(realIx)]);
            CovZP = cov([xDepth{rh,jj,ii}(realIx)],[xRhoRef{rh,jj,ii}(realIx)]);
            HorizonDepthVar{rh,jj,ii} = CovZP(1,1);
            RhoRefVar{rh,jj,ii} = CovZP(2,2);
            CovDepthRho{rh,jj,ii} = CovZP(1,2);
            
        end
        
        % Perform Dix Inversion of RMS velocities
        parfor (jj = looper, nWorkers)
%         for jj = looper
            if rh > 1
                % Real Solution Constraint
                realIxI = find(real([xToRef{rh-1,jj,ii}]));
                realIxJ = find(real([xToRef{rh,jj,ii}]));
                
                [xVint{rh,jj,ii}, xHint{rh,jj,ii}] ...
                    = DixHVA([xVref{rh-1,jj,ii}(realIxI)],[xVref{rh,jj,ii}(realIxJ)],...
                    [xToRef{rh-1,jj,ii}(realIxI)],[xToRef{rh,jj,ii}(realIxJ)]);
                xRhoInt{rh,jj,ii} = DryCrim(xVint{rh,jj,ii});
                VoInt{rh,jj,ii} = mean([xVint{rh,jj,ii}]);
                VoIntVar{rh,jj,ii} = var([xVint{rh,jj,ii}]);
                HoInt{rh,jj,ii} = mean([xHint{rh,jj,ii}]);
                HoIntVar{rh,jj,ii} = var([xHint{rh,jj,ii}]);
                RhoInt{rh,jj,ii} = mean([xRhoInt{rh,jj,ii}]);
                RhoIntVar{rh,jj,ii} = var([xRhoInt{rh,jj,ii}]);
                CovHP = cov([xHint{rh,jj,ii}],[xRhoInt{rh,jj,ii}]);
                CovThicknessRho{rh,jj,ii} = CovHP(1,2);
                
            else
                realIx = find(real([xToRef{rh,jj,ii}]));
                xVint{rh,jj,ii} = xVref{rh,jj,ii}(realIx);
                xHint{rh,jj,ii} = xDepth{rh,jj,ii}(realIx);
                xRhoInt{rh,jj,ii} = xRhoRef{rh,jj,ii}(realIx);
                HoInt{rh,jj,ii} = HorizonDepth{rh,jj,ii};
                HoIntVar{rh,jj,ii} = HorizonDepthVar{rh,jj,ii};
                VoInt{rh,jj,ii} = VoRef{rh,jj,ii};
                VoIntVar{rh,jj,ii} = VoVarRef{rh,jj,ii};
                RhoInt{rh,jj,ii} = RhoRef{rh,jj,ii};
                RhoIntVar{rh,jj,ii} = RhoRefVar{rh,jj,ii};
                CovHP = cov([xHint{rh,jj,ii} ],[xRhoInt{rh,jj,ii}]);
                CovThicknessRho{rh,jj,ii} = CovHP(1,2);
            end
            
        end

                    
        % Concatenate Reflection Travel Time
        ReflectionTo{rh,ii} = [toRef{rh,:,ii}];
        ReflectionToVar{rh,ii} = [toVarRef{rh,:,ii}];
        % Concatenate Reflection Velocity
        ReflectionVelocity{rh,ii} = [VoRef{rh,:,ii}];
        ReflectionVelocityVar{rh,ii} = [VoVarRef{rh,:,ii}];
        IntervalVelocity{rh,ii} = [VoInt{rh,:,ii}];
        IntervalVelocityVar{rh,ii} = [VoIntVar{rh,:,ii}];
        % Concatenate Reflector Depth
        ReflectionDepth{rh,ii} = [HorizonDepth{rh,:,ii}];
        ReflectionDepthVar{rh,ii} = [HorizonDepthVar{rh,:,ii}];
        %Concatenate Layer Thicnkess
        IntervalThickness{rh,ii} = [HoInt{rh,:,ii}];
        IntervalThicknessVar{rh,ii} = [HoIntVar{rh,:,ii}];
        % Concatenate Reflection Density
        ReflectionDensity{rh,ii} = [RhoRef{rh,:,ii}];
        ReflectionDensityVar{rh,ii} = [RhoRefVar{rh,:,ii}];
        % Concatenate Interval Density
        IntervalDensity{rh,ii} = [RhoInt{rh,:,ii}];
        IntervalDensityVar{rh,ii} = [RhoIntVar{rh,:,ii}];
        % Concatenate Covariance
        CovarianceDepthDensity{rh,ii} = [CovDepthRho{rh,:,ii}];
        CovarianceThicknessDensity{rh,ii} = [CovThicknessRho{rh,:,ii}];
        
        % Estimte SWE
        % Initiate CallBack to GreenTracsFirnCore
        if isGreenTracsFirnCore
            isCallbackHVA = 1;
            % Assumes Firn Stratigraphy to be Isochronous
            % Report the Age of the Reflection Horizon 
            GreenTracsFirnCore
        end
        % Estimate Annual Accumulation of ageInterval Annuals
        SWE{rh,ii} = ReflectionDepth{rh,ii}.*ReflectionDensity{rh,ii}./ageInterval;
        SWEint{rh,ii} = IntervalThickness{rh,ii}.*IntervalDensity{rh,ii}./ageInterval;
        % Error Propagation Equation
        SWEvar{rh,ii} = ReflectionDepthVar{rh,ii}.*ReflectionDensity{rh,ii}.^2 ...
            + ReflectionDensityVar{rh,ii}.*ReflectionDepth{rh,ii}.^2 ...
            + 2.*ReflectionDepth{rh,ii}.*ReflectionDensity{rh,ii}.*CovarianceDepthDensity{rh,ii};
        SWEintVar{rh,ii} = IntervalThicknessVar{rh,ii}.*IntervalDensity{rh,ii}.^2 ...
            + IntervalDensityVar{rh,ii}.*IntervalThickness{rh,ii}.^2 ...
            + 2.*IntervalThickness{rh,ii}.*IntervalDensity{rh,ii}.*CovarianceThicknessDensity{rh,ii};
        
        % Smooth Esitmates
        smoothR = 251;
        if isReduceData
            if mod(smoothR,2) == 0
            smoothR = ceil(smoothR./rmNtrc);
            else
                smoothR = smoothR+1;
            end
        end
        TWT{rh,ii} = nonParametricSmooth( 1:length(ReflectionTo{rh,ii}),...
            ReflectionTo{rh,ii},1:length(ReflectionTo{rh,ii}),smoothR);
        TWTvar{rh,ii} = nonParametricSmooth( 1:length(ReflectionToVar{rh,ii}),...
            ReflectionToVar{rh,ii},1:length(ReflectionToVar{rh,ii}),smoothR);
        SnowWaterEqv{rh,ii} = nonParametricSmooth( 1:length(SWE{rh,ii}),SWE{rh,ii},...
            1:length(SWE{rh,ii}),smoothR);
        SnowWaterEqvVar{rh,ii} = nonParametricSmooth( 1:length(SWEvar{rh,ii}),...
            SWEvar{rh,ii},1:length(SWEvar{rh,ii}),smoothR);
        Density{rh,ii} = nonParametricSmooth( 1:length(ReflectionDensity{rh,ii}),...
            ReflectionDensity{rh,ii},1:length(ReflectionDensity{rh,ii}),smoothR);
        DensityVar{rh,ii} = nonParametricSmooth( 1:length(ReflectionDensityVar{rh,ii}),...
            ReflectionDensityVar{rh,ii},1:length(ReflectionDensityVar{rh,ii}),smoothR);
        Depth{rh,ii} = nonParametricSmooth( 1:length(ReflectionDepth{rh,ii}),...
            ReflectionDepth{rh,ii},1:length(ReflectionDepth{rh,ii}),smoothR);
        DepthVar{rh,ii} = nonParametricSmooth( 1:length(ReflectionDepthVar{rh,ii}),...
            ReflectionDepthVar{rh,ii},1:length(ReflectionDepthVar{rh,ii}),smoothR);
        % Smooth Inteval Estimates
        LayerSnowWaterEqv{rh,ii} = nonParametricSmooth( 1:length(SWEint{rh,ii}),SWEint{rh,ii},...
            1:length(SWEint{rh,ii}),smoothR);
        LayerSnowWaterEqvVar{rh,ii} = nonParametricSmooth( 1:length(SWEintVar{rh,ii}),...
            SWEintVar{rh,ii},1:length(SWEintVar{rh,ii}),smoothR);
        LayerDensity{rh,ii} = nonParametricSmooth( 1:length(IntervalDensity{rh,ii}),...
            IntervalDensity{rh,ii},1:length(IntervalDensity{rh,ii}),smoothR);
        LayerDensityVar{rh,ii} = nonParametricSmooth( 1:length(IntervalDensityVar{rh,ii}),...
            IntervalDensityVar{rh,ii},1:length(IntervalDensityVar{rh,ii}),smoothR);
        LayerThickness{rh,ii} = nonParametricSmooth( 1:length(IntervalThickness{rh,ii}),...
            IntervalThickness{rh,ii},1:length(IntervalThickness{rh,ii}),smoothR);
        LayerThicknessVar{rh,ii} = nonParametricSmooth( 1:length(IntervalThicknessVar{rh,ii}),...
            IntervalThicknessVar{rh,ii},1:length(IntervalThicknessVar{rh,ii}),smoothR);
        end
    end
