%% GreenTracsFirnCore 
% Uses laboratory measurements of firn age-depth scales and accumulation
% to improve the accuracy of the MxHL model.
%
% Citations Needed

% Loop For Additional Cores Coming Soon
%% Receive Callback from Horizon Velocity Analysis
try isCallbackHVA;
catch
    isCallbackHVA = 0;
end
    
if isCallbackHVA
    [~,horizonIx] = min(abs(GTCdepthAge(:,1)...
        - mean(ReflectionDepth{rh,ii}(iceCoreIx))));
    ageInterval = GTCdepthAge(1,2)-GTCdepthAge(horizonIx,2);
else
%% Load 2-column format .txt [Lat, Lon; depth, year]
GTCza = load(depthAgeFilename);
GTCzp = load(depthDensityFilename);

%% Extract GPS Location of Firn Core
iceCoreLL = GTCza(1,:);
iceCoreXY = ll2psn(iceCoreLL(1),iceCoreLL(2));
if norm(GTCzp(1,:) - iceCoreLL) < 1e-8
GTCza(1,:) = [];
GTCzp(1,:) = [];
elseif norm(GTCzp(1,:) - iceCoreXY) < 1e-8
GTCza(1,:) = [];
GTCzp(1,:) = [];
else
    error('Age and density data are not coloacted, or are in the wrong format. Check the supplemental data!');
end

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

%% Interpolate Depth-Age Scale to Model Resolution
% Append Snow Surface Date to Age Record
GTCza = [0,str2num(Year{1})+dayofyear/365;GTCza];
GTCaz = fliplr(GTCza);


% Array of Model Depths
maxDepth = max(GTCza(:,1));
maxAge = min(GTCza(:,2));
modDepth = 30.001;
nDepth = 301;
tmpZ = (linspace(0.001,modDepth,nDepth))';

% Determine Model Depths Greater than Core
[~, xstrapIx] = min(abs(tmpZ-maxDepth));

% Interpolate Age-Depth Scale
% depthAge : Depth is the linear Axis!
GTCdepthAge = pchip(GTCza(:,1),GTCza(:,2),tmpZ(1:xstrapIx-1));

% Extrapolate Age-Depth with Linear Model 
% Estimating slope from previous 1 year of accumulation allows
% extrapolation to be consistent with the core and HL1980 model.
[~, modIx] = min(abs(GTCdepthAge-maxAge - 1));
modIx = modIx:xstrapIx-1;

% L2 Norm Regression
G = [ones(length(modIx),1),tmpZ(modIx)];
d = [GTCdepthAge(modIx)];
m = G\d;

% Extrapolate Ages to 30 m
GTCdepthAge = [GTCdepthAge;m(2).*tmpZ(xstrapIx:end)+m(1)];
GTCdepthAge = [tmpZ,GTCdepthAge];
clear('modIx','G','m','d','xstrapIx','tmpZ','nDepth','GTCza',...
    'maxAge','maxDepth','modDepth','depthAgeFilename')

%% Load GreenTrACS Accumulation Record
Accumulation = load('Accumulation.txt');

n = 10;                         % Years to be sampled
N = 1250;                       % N Monte Carlo simulations
nCore = length(unique(coreNo)); % Number of Cores in Analysis
Accu = zeros(N,nCore);          % Allocation

% Monte Carlo Simulation for Average Annual Accumulation
for ii = 1:nCore
    nanIx = isnan(Accumulation(:,coreNo(ii)));
    for jj = 1:N
        Accu(jj,ii) = mean(datasample(Accumulation(~nanIx,coreNo(ii)),n,'replace',false));
    end
end
% Average Value from the Central Limit
coreAccumulation = mean(Accu);
coreAccumulation(2,:) = std(Accu);
% Clear the Memory
clear('Accumulation','Accu','n','N','nanIx')
end


