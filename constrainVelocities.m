function [ NewV, NewTo, NewZ, Vint,Hint ] = constrainVelocities( ReflectionFBpick, AirTo, xToDir,xToRef, xVdir, xVref, xDepthDir, xDepth, offsetArray)
% Inputs:   
%           Vlmo - xVdir
%           Vrms - xVref
%           Zlmo - xDepthDir
%           Znmo - xDepth
%           Tnmo - ReflectionFBpick
%           ToAir- toDir{1}
%           offsetArray
%
% Output:   newVrms
%           newZnmo
%           newTo
%           Vint?
%           Hint?

% Allocate BootStrapped Distributions
xVlmo = cell2mat(xVdir');
xVnmo = cell2mat(xVref');
xZlmo = cell2mat(xDepthDir');
xZnmo = cell2mat(xDepth');
xTolmo = abs(cell2mat(xToDir'));
xTonmo = cell2mat(xToRef');

% Catch Bad Sampling
isBadSample = 1;
ii = 1;
% Run Monte Carlo Simulation
while ii <= 5;%250
    isBadSample = 1;
    while isBadSample
        % Randomply Sample Bootstrapped Distributions
%         SampleIx = 1:size(xVnmo,1);
        randIx = 1:size(xVnmo,2);
%         randIx = randsample(SampleIx,size(xVnmo,2));
%         [Vlmo,lmoIx] = datasample(xVlmo,1);
%         for jj = 1:length(randIx)
%             Vnmo(jj) = xVnmo(randIx(jj),jj);
%             Znmo(jj) = xZnmo(randIx(jj),jj);
%             Tonmo(jj) = xTonmo(randIx(jj),jj);
%         end

        for jj = 1:length(randIx)
            Vnmo(jj) = mean(xVnmo(:,jj));
            Znmo(jj) = mean(xZnmo(:,jj));
            Tonmo(jj) = mean(xTonmo(:,jj));
        end
        Vlmo = mean(xVlmo);Vlmo = Vlmo(:);
        Zlmo = mean(xZlmo);Zlmo = Zlmo(:);
        Tolmo = mean(xTolmo);Tolmo = Tolmo(:);
        Vnmo = Vnmo(:);
        Znmo = Znmo(:);
        Tonmo = Tonmo(:);
        % [Vnmo,nmoIx] = datasample(xVnmo,1); Vnmo = Vnmo(:);
%         Zlmo = xZlmo(lmoIx,:);
%         Znmo = xZnmo(nmoIx,:); Znmo = Znmo(:);
%         Tolmo = xTolmo(lmoIx,:);
        To = [Tolmo;Tonmo];
%         Toint = [Tolmo;diff(To./2)];
        ToInt = diff(To./2);
        Zo = [Zlmo;Znmo];
        % Constrain Depth Gradient
        ZoSort = sort(Zo);
        if any(Zo ~= ZoSort)
            isBadSample = 1;
            ii = ii - 1;
            break
        else isBadSample = 0;
        end
        
        
        
        % Initialize Loop State
        rmvHorizon = 0;
        % Testable Gradient Sigma Tolerance
        tol = 3;%1:.5:100;
        % Tolerance Iterator
        % tt = 1;
        
        % User May Sacrifice Horizons to Increase Smoothness Tolerance
        % If Smoothness is Sacrificed (Default) the Algorithm Applies the Minumum
        % Tolerance which accpets all input data
        %   - This Approach Comes at Computational Cost
        
        % Loop Through Tolerances if Horizons are not Sacrificed
%         while ~isempty(rmvHorizon);
        % Need to Input Data From SWEDISH here
        do = [Vlmo;Vnmo];
        d = do(2:end);
        H = [Zlmo;diff([Zlmo;Znmo])];
        
        % Initialize Vector of Horizon Numbers - 1 is the Surface
        HorizonIx = 1:length(H);
        % clear rmvHorizon
        rmvHorizon = [];
        
        % Impose Herron & Langway Model for Invertable/Stable Vrms Gradient
        % Depths for Model Estimate are Solved from Vrms Profile
%         Zhl = [Zlmo;Znmo];
        Zhl = linspace(0.001,30.001,301);Zhl = Zhl(:);
        Hhl = [0.001;diff(Zhl)];
        % Assume Mean Annual Temperature at -20C
        annualT = -20;
        % Compute Surface Density from Direct Wave Velocity
        surfDensity = DryCrim(do(1));
        % Assume Annual Accumulation at 0.5 m.w.e.
        annualA = 0.5;
        % Compute HL Model Density
        densityHL = herronLangway( Zhl, annualT, surfDensity, annualA);
        % Compute Cumulative Average Density
        cumW = zeros(length(densityHL),length(densityHL));
        for jj = 1:length(densityHL)
            cumW(jj,1:jj) = Hhl(1:jj)./Zhl(jj);
        end
        densityAvgHL = cumW*densityHL;
%         densityAvgHL = cumsum(densityHL)./(1:length(densityHL))';
        % Estimate Vrms from HL model ~ to Average Density
        Vhl = DryCrimVRMS(densityAvgHL);
%         % Calculate Zero-offset Times from HL NMO Velocity
%         bumhl = zeros(size(Vhl));
%         for jj = 1:length(Vhl)-1
%             % Cast this into a Matrix
%             xs = (offsetArray.^2.*(Vhl(jj+1).^(-2)))';
%             x1 = ones(length(xs),1);
%             d2 = (ReflectionFBpick(:,jj)-AirTo).^2;
%             Ghlto = [x1,xs];
%             mto = Ghlto\d2;
%             bumhl(jj+1) = sqrt(mto(1));
%             %                         toNMO(forceIx(jj)-1) = sqrt(mto(1));
%         end
% %         for jj = 1:length(randIx)
% %             toHL(jj) = mean(sqrt(abs((ReflectionFBpick(:,jj)-AirTo).^2 - (offsetArray.^(2)./(Vhl(jj+1).^2))')));
% %         end
% newZhl = bumhl(2:end).*Vhl(2:end)./2;
% newZhl = [Zhl(1);newZhl];
% errZhl = Zhl-newZhl;
% newH = [newZhl(1);diff(newZhl)];
%         % Compute HL Model Density
%         newDensityHL = herronLangway( newZhl, annualT, surfDensity, annualA);
%         % Compute Cumulative Average Density
%         cumW = zeros(length(newDensityHL),length(newDensityHL));
%         for jj = 1:length(newDensityHL)
%             cumW(jj,1:jj) = newH(1:jj)./newZhl(jj);
%         end
%                 newDensityAvgHL = cumW*densityHL;
% %         densityAvgHL = cumsum(densityHL)./(1:length(densityHL))';
%         % Estimate Vrms from HL model ~ to Average Density
%         newVhl = DryCrimVRMS(newDensityAvgHL);
        
        toHL = Zhl./Vhl;%toHL(1) = toHL(1);
        toHL = toHL(:); %toHL = [Tolmo;toHL];ToIntHL = diff(toHL./2);
        % Find HL Index Nearest the Radar Estimation
        [~,HLix] = min(pdist2(To./2,toHL),[],2);
        ToIntHL = diff(toHL(HLix));
        %one-way Time or two-way time? methinks one-way-time
        % Dix solution of HL model
%         GhlInt = tril(([toHL(1);ToIntHL]*ones(1,length([toHL(1);ToIntHL])))');
%         VrmsHL = (Vhl.^2).*toHL;
%         mVintHL = GhlInt\VrmsHL;
%         VintHL = sqrt(mVintHL);
%         HintHL = VintHL(:,ii).*[toHL(1);ToIntHL];
%        
        
        % Constrained RMS Gradient Inversion
        % Model G is the Causal Integral Operator
%         G = tril(ones(length(d)));
        % Estimate of the Average Velocity Gradient
%         x = gradient(d,Toint);
% x = gradient(d);
        % Gradients not right here.
        x = diff(do)./(ToInt);
%         x = diff(do)./(To(2:end)./2);
% 12:00 Today
% x = do./To./2;
% xHL = gradient(Vhl(2:end));
        xHL = diff(Vhl(HLix))./(ToIntHL);
% xHL = diff(Vhl)./(toHL);
%         xHL = diff(Vhl)./(toHL(2:end)./2);
ToHL = toHL;
        clear toHL
%         xHL = gradient(Vhl);
%         x = G\d;    % Data Estimate
%         xHL = G\Vhl;% HL Estimate
        % Calculate the Average and Variance of the HL Model to Guide Vrms update
        HLgradient = mean(xHL);
        % HLgradVar = tol(tt).*std(xHL(2:end));
        HLgradVar = tol.*std(xHL);
        
        % Stable Vrms Gradient Follows the Herron-Langway Gradient within Tolerance
        % Calcuate the Residual Velocity Gradient
        residual = (x-(xHL));
%         h = find(x > HLgradient + HLgradVar | x < HLgradient - HLgradVar,1);
        h = find(x > xHL + HLgradVar | x < xHL - HLgradVar,1);

        iter = 0;
        while x(h) > xHL(h) + HLgradVar || x(h) <= xHL(h) - HLgradVar
%         while x(h) > HLgradient + HLgradVar || x(h) <= HLgradient - HLgradVar
            
            Nudge = abs(x(h)+0.000005).*(ToInt(h)./5);
%             Nudge = abs(HLgradient+0.000005).*(1./2);

            % If gradient is Positive Nudge is Negative
%             if x(h) > HLgradient + HLgradVar
            if x(h) > xHL(h) + HLgradVar
                Nudge = -Nudge;
            end
            if abs(Nudge) < 0.0005
                Nudge = sign(Nudge)*0.0015;
            end

            d(h) = d(h)+(Nudge);
            do(h+1) = d(h);
            
        % Calculate Zero-offset Times from new NMO Velocity
        toNMO = zeros(length(randIx),1);
        for jj = 1:length(randIx)
            % Cast this into a Matrix
            xs = (offsetArray.^2.*(d(jj).^(-2)))';
            x1 = ones(length(xs),1);
            d2 = (ReflectionFBpick(:,jj)-AirTo).^2;
            Gto = [x1,xs];
            mto = Gto\d2;
            toNMO(jj) = sqrt(mto(1));
        end
%         for jj = 1:length(randIx)
%             toNMO(jj) = mean(sqrt((ReflectionFBpick(:,jj)-AirTo).^2 - (offsetArray.^(2)./(d(jj).^2))'));
%         end
%         toNMO = toNMO(:); 
        To = [Tolmo;toNMO];ToInt = diff(To./2);
        % ReCalculate Updated Gradient
            x = diff(do)./(ToInt);
                                % Find HL Index Nearest the Radar Estimation
                    [~,HLix] = min(pdist2(To./2,ToHL),[],2);
                    ToIntHL = diff(ToHL(HLix));
            xHL = diff(Vhl(HLix))./(ToIntHL);
%                     % Find HL Index Nearest the Radar Estimation
%         [~,HLix] = min(pdist2(To./2,toHL),[],2);
%         ToIntHL = diff(toHL(HLix));
%             % ReCharacterize Horizon h

%             newh = find(find(x > HLgradient + HLgradVar | x < HLgradient - HLgradVar)>4,1);
%             newh = find(x > HLgradient + HLgradVar | x < HLgradient - HLgradVar,1);
            newh = find(x > xHL + HLgradVar | x < xHL - HLgradVar,1);
            if newh == h
                iter = iter + 1;
                if iter>50
                    % Find HL Index Nearest the Radar Estimation
                    [~,HLix] = min(pdist2(To(h)./2,ToHL),[],2);
                    d(h) = Vhl(HLix);
                    do(h+1) = d(h); 
                    iter = 0;
                    x(h) = xHL(h);
                    h = find(x > xHL + HLgradVar | x < xHL - HLgradVar,1);
                end
            else h = newh;
            end
            if isempty(h)
                % Find HL Index Nearest the Radar Estimation
                [~,HLix] = min(pdist2(To./2,ToHL),[],2);
                newVhl = Vhl(HLix);
                % Check Residuals
                residualV = (do-newVhl);
                % Force Significant Residual to HL Model
                forceIx = find(abs(residualV) > 0.0005);
                if forceIx(1)==1
                    forceIx(1) = [];
                end
                % Perturbated HL Model Solution
                do(forceIx) = newVhl(forceIx);% + 20.*Vhl(forceIx).*x(forceIx-1);
                %                     do(forceIx) = Vhl(forceIx);
                % Re-Calculate Zero-offset Times from new NMO Velocity
                %                     toNMO = zeros(length(randIx),1);
                for jj = 1:length(forceIx)
                    % Cast this into a Matrix
                    xs = (offsetArray.^2.*(do(forceIx(jj)).^(-2)))';
                    x1 = ones(length(xs),1);
                    d2 = (ReflectionFBpick(:,forceIx(jj)-1)-AirTo).^2;
                    Gto = [x1,xs];
                    mto = Gto\d2;
                    %                         bum(forceIx(jj)-1) = sqrt(mto(1));
                    toNMO(forceIx(jj)-1) = sqrt(mto(1));
                end
                break
                %                 end
            end
        end
        % Interval Velocity Inversion
        NewV(:,ii) = do;
        NewTo(:,ii) = [Tolmo;toNMO];
        NewZ(:,ii) = do.*NewTo(:,ii)./2;
        NewT1 = [NewTo(1,ii);NewTo(2:end,ii)./2];
        NewIntT = [NewT1(1);diff(NewT1)];
        Gint = tril((NewIntT*ones(1,length(NewIntT)))');
        Vrms2 = (do.^2).*NewT1;
        mVint = Gint\Vrms2;
        Vint(:,ii) = sqrt(mVint);
        Hint(:,ii) = Vint(:,ii).*NewT1;
        
    end
%     ii = ii + 1;
ii = 1000;
end

% Physically Constrain Interval Velocities
[~,imagJx] = find(imag(Vint));
[~,falseJx] = find(abs(Vint)>0.25 | abs(Vint)<0.1689);
RemoveJx = unique([imagJx;falseJx]);
Vint(:,RemoveJx) = [];
Hint(:,RemoveJx) = [];
NewV(:,RemoveJx) = [];
NewTo(:,RemoveJx) = [];
NewZ(:,RemoveJx) = [];
% Export Ensemble Average
Vint = mean(Vint,2);
VintVar = var(Vint,0,2);
Hint = mean(Hint,2);
HintVar = var(Hint,0,2);
NewV = mean(NewV,2);
NewVvar = var(NewV,0,2);
NewTo = mean(NewTo,2);
NewToVar = var(NewTo,0,2);
NewZ = mean(NewZ,2);
NewZvar = var(NewZ,0,2);






end