%% Read Sensors and Software Data
% This Script Reads the Binary Multi-channel GPR Data and Performs the
% pre-processing and signal processing routines. 
% GPS Locations are incorporated in NSIDC Sea Ice Polarstereographic North
%
% The Array Geometry is Installed
% The Acquisition Method is Decided
% The Near Channel (7) is Killed
% The Time Window can be trimmed
% Data is de-multiplex (grouped into channels) and coarse trace shifts
% are applied as the first cut at correcting for digital time-sampling 
% errors.
% Inside the function wrapper processCommonOffset.m de-WOW filtering and
% trace stacking are applied. There are Various filter parameters decided
% in this function.

    Rad = cell(1,nFiles);
    trhd = cell(1,nFiles);
%     GPS = cell(1,nFiles);
    Year = cell(1,nFiles);
    TimeAxis = cell(1,nFiles);
    
    for ii = 1 : nFiles
        tic
        %------------------------------------------------------------------
        % Multiplexed Channel Record
        filename = fileNames(lineNo(ii)+1).name;
        filepath = fullfile(dataDir,filename);
        % Read netCDF data file
        disp(' ')
        fprintf('Reading .nc File \n')
        tic
        ncRad = ncread(filepath,'DATA');
        trhd{ii} = ncRad(1:25,:);
        Rad{ii} = ncRad(26:end,:);
        clear('ncRad');
        f0 = (trhd{ii}(9,1)); % [MHz]
        dt = trhd{ii}(7,1)/1000; % [ns]
        dx = diff(trhd{ii}(2,1:2)); % [m]
        % Need to Automate Offset Array from .nc File
%         offsetArray = unique(trhd{ii}(22,1:100));

        % Configure GPS
        if isLoadGPS
            disp(' ')
            fprintf('Configuring GPS \n')
            tic
            % GPS = [X,Y,Z,Distance,Slope,Speed,Heading,Tailing];
            GPSixEdges = GeoLocation.GPSixEdges(:,ii);
            GPSix = GPSixEdges(1):GPSixEdges(2);
            Year{ii} = GeoLocation.Date{GPSix(1)}(1:4);
            tmpDate = datetime(GeoLocation.Date{GPSix(1)},'InputFormat','yyyy/MM/dd');
            dayofyear = day(tmpDate,'dayofyear');
            GPS = [GeoLocation.X(GPSix),GeoLocation.Y(GPSix),...
                GeoLocation.Z_EGM08(GPSix)];
            dX = GeoLocation.dX(GPSix);
            dY = GeoLocation.dY(GPSix);
            dZ = GeoLocation.dZ(GPSix);
            
            % Compute Distance
            dS = sqrt(dX.^2+dY.^2+dZ.^2);
            if ii>1
                Distance = endDist+cumsum(dS);
            else
                Distance = cumsum(dS); %[m]
            end
            endDist = Distance(end);
%             DistanceKm = Distance./1000; % [km]

            % Slope [Degrees]
            Slope = GeoLocation.Slope(GPSix);
            
            % Speed [m/s]
            Speed = GeoLocation.Speed(GPSix);
            
            % Heading [Degrees]
            Heading = GeoLocation.Heading(GPSix);
            Tailing = mod(Heading+180,360);
            GPS = [GPS,Distance,Slope,Speed,Heading,Tailing];
            clear('GPSix','GPSixEdges','dS','dX','dY','dZ',...
                'Distance','Slope','Speed','Heading','Tailing');
        end
 
        % Nominal Frequency GHz
        f0GHz = f0/1000;
        % No. Traces in Multiplexed Data
        [~, multiplexNtrcs] = size(Rad{ii});
        
        % Install Transmitter and Receiver Sequencing and Geometry
        % 500 MHz Offsets
        if f0 == 500
            if nChan == 9
                % 9 Channel Offsets
                txGeo = [0, 0, 0, 1.33, 1.33, 1.33, 2.67, 2.67, 2.67]; % Tx Sequence & Absolute Position
                rxGeo = [4, 8, 12, 4, 8, 12, 4, 8, 12]; % Rx Sequence & Absolute Position
            elseif nChan == 6 && ~isKill
                % 6 Channel Offsets
                txGeo = [1.33, 1.33, 1.33, 2.67, 2.67, 2.67]; % Tx Sequence & Absolute Position
                rxGeo = [4, 8, 12, 4, 8, 12]; % Rx Sequence & Absolute Position
            end
            % 1 GHz Offsets
        elseif f0 == 1000
            txGeo = [0, 0, 0, 0.67, 0.67, 0.67, 1.33, 1.33, 1.33]; % Tx Sequence & Absolute Position
            rxGeo = [2, 4, 6, 2, 4, 6, 2, 4, 6]; % Rx Sequence & Absolute Position
        end
            
        offsetArray = rxGeo - txGeo;            % Offset Array for CMP Gather
        midPointArray = mean([txGeo;rxGeo],1);  % Midpoint Locations Relative to Tx 1
        trhd{ii}(21,:) = midPointArray([trhd{ii}(23,:)]); % Append Midpoint Locations
        nChan = length(offsetArray);            % Refresh nChan before kill
        farOffset = max(offsetArray);           % Determine Far Offset
        farChan = find(offsetArray == farOffset);
        nearOffset = min(offsetArray);          % Determine Near Offset
        nearChan = find(offsetArray == nearOffset);
        
        % Determine Data Acquisition  Method
        if all(trhd{ii}(5,:) == 0)
            isFreeRun = 1;
            isOdometer = 0;
        else
            isOdometer = 1;
            isFreeRun = 0;        
        end
        
        % Remove Static Traces if is Free Run Acquisition
        if isFreeRun
            % Process Near Offset Channel
            disp(' ')
            fprintf('Begin Static Trace Removal \n')
            tic
            nearRad = processTraceRemoval... % Processes Near Offset Data
                (Rad{ii}(:,nearChan:nChan:end), f0, dt );
            % Remove Duplicate Static Traces
            dupIx = removeStaticTrace( nearRad, multiplexNtrcs, nearChan, nChan );
            clear('nearRad');
            
            trhd{ii}(:,dupIx) = [];         % Remove Static Trace Headers from Multiplexed Record
            trhd{ii}(1,:) = 1:length(trhd{ii}); % Configure Trace Indicies
            trhd{ii}(2,:) = [0:dx:(length(trhd{ii})-1)*dx]; % Configure Distance
            Rad{ii}(:,dupIx) = [];      % Remove Static Traces from Multiplexed Data
            xArray = trhd{ii}(2,:);         % Define Configured Distance xArray
            
            fprintf('Static Trace Removal Done \n')
            toc
            display(' ')
            
        end
        
        % Interpolate GPS
        % GPS = [X,Y,Z,Distance,Slope,Speed,Heading,Tailing];
        nGPS = size(GPS,1);
        nTrcs = size(Rad{ii},2);
        xq = linspace(1,nGPS,nTrcs);
        % Piecewise Cubic Hermite Interpolating Polynomials
        Si = pchip(1:nGPS,GPS(:,4),xq);
        Xi = pchip(GPS(:,4),GPS(:,1),Si);
        Yi = pchip(GPS(:,4),GPS(:,2),Si);
        Zi = pchip(GPS(:,4),GPS(:,3),Si);
        Sxi= pchip(GPS(:,4),GPS(:,5),Si);
        Vi = pchip(GPS(:,4),GPS(:,6),Si);
        Hi = pchip(GPS(:,4),GPS(:,7),Si);
        Ti = pchip(GPS(:,4),GPS(:,8),Si);
        % Append GPS to Trace Header
        trhd{ii}(13:20,:) = [Xi;Yi;Zi;Si;Sxi;Vi;Hi;Ti];
        clear('GPS','nGPS','nTrcs','Si','Xi','Yi','Zi','Sxi','Vi','Hi','Ti','xq');
        
        % Remove Skip Traces if Wheel Odometer Acquisition
        if isOdometer
        dupIx = find(~(diff(trhd{ii}(23,:)))); % Find Skipped Traces
        
        for jj = dupIx
            trhd{ii}(1,jj:end) = trhd{ii}(1,jj:end) - 1; % ReConfigure Trace Indicies
        end
        
        trhd{ii}(:,dupIx) = []; % Remove Skipped Traces from Trace Header
        Rad{ii}(:,dupIx) = []; % Remove Skipped Traces from Multiplexed Data
        xArray = trhd{ii}(2,:); % Define ReConfigured Distance as xArray
        end
        % Gain Far Channels for Cable Attenuation
        if f0 == 500
            gainIx = find(trhd{ii}(23,:) == 3 | trhd{ii}(23,:) == 6 | trhd{ii}(23,:) == 9);
            Rad{ii}(:,gainIx) = Rad{ii}(:,gainIx).*1.7;
        end
                
        % Kill Unwanted Channels
        if isKill
            if ii == 1
                killChan = [7];         % Kill Channel Number
                liveChan(killChan) = [];% Array of Live Channels
                % Update Killed Channel Shifts
                chanShift(killChan) = [];
            end
            killIx = find(ismember(trhd{ii}(23,:),killChan));
            % Kill Traces
            Rad{ii}(:,killIx) = [];  
            % Remove Killed Trace Headers
            trhd{ii}(:,killIx) = [];        
            % Configure Trace Indicies
            trhd{ii}(1,:) = 1:length(trhd{ii}); 
            % Define Configured Distance xArray
            xArray = trhd{ii}(2,:);        
            % Remove Killed Offsets
            offsetArray(killChan) = []; 
            % Number of Live Channels
            nChan = nChan - length(killChan);
            % Looping Array of Record Channels
            chan =  1:nChan;
            % Determine Far Offset After Trace Kill
            farOffset = max(offsetArray);           
            farChan = find(offsetArray == farOffset);
            % Determine Near Offset After Trace Kill
            nearOffset = min(offsetArray);          
            nearChan = find(offsetArray == nearOffset);
        end
        clear('dupIx','gainIx','killIx');

        % Pad Data with Instrument Zero
        if strcmp(Year{ii},'2016')
            padding = 50;
        end
        if strcmp(Year{ii},'2017')
            padding = 0;
        end
        instrumentPad = zeros(padding,size(Rad{ii},2));
        if padding ~= 0
            for jj = 1:size(Rad{ii},2)
                instrumentZero = Rad{ii}(1,jj);
                instrumentPad(:,jj) = ones(padding,1).*instrumentZero;
            end
            Rad{ii} = [instrumentPad;Rad{ii}];
        end  
        
        % Trim Time Window
        if isTrimTWT
            reSample = 575 + padding;   % Number of Wanted Samples
            Rad{ii} = Rad{ii}(1:reSample,:);
        end
        
        % Allocation Here
        if ii == 1
            Radar = cell(nChan,nFiles); traceIx = cell(nChan,nFiles);
%             Array = cell(nChan,nFiles);
        end
        
        for jj = chan
            % DeMux Sequential Data
            % GPS DeadReckoning Within Demux
                % Calculate Truer Antenna Positions using Array Geometry
            [Radar{jj,ii},trhd{ii},traceIx{jj,ii},~] = DeMux(Rad{ii},trhd{ii},liveChan(jj),isLoadGPS);
        end
        
        multiplexNtrcs = size(trhd{ii},2);% Length of Multiplex
        traceMod = mod(multiplexNtrcs,nChan);% Count Unbinned Traces
        trhd{ii}(:,(multiplexNtrcs - traceMod)+1 :multiplexNtrcs ) = []; % Remove Unbinned Trace Headers

        % Average Positions for Bin Centers (The GPS Location and Distance)
        if isLoadGPS
            for jj = 1:nChan:size(trhd{ii},2)
                % Store Bin Centers [m]
                trhd{ii}(10:12,jj:jj+nChan-1) = mean(trhd{ii}(13:15,jj:jj+nChan-1),2)*ones(1,nChan);
                % Overwrite Distance with Bin Center Position [m]
                tmp = mean(trhd{ii}(16,jj:jj+nChan-1),2)*ones(1,nChan);
                trhd{ii}(2,jj:jj+nChan-1) = tmp;
                % Overwrite Tailing with Average Bin Center Heading
                trhd{ii}(20,jj:jj+nChan-1) = mean(trhd{ii}(19,jj:jj+nChan-1),2)*ones(1,nChan);
            end
            % Smooth Heading
%             trhd{ii}(20,:) = nonParametricSmooth(1:length(trhd{ii}),trhd{ii}(20,:),1:length(trhd{ii}),251.*nChan);
            % Zero Starting Distance
            tmp = trhd{1}(2,1);
            trhd{ii}(2,:) = trhd{ii}(2,:) - tmp;
            
        end
        clear('tmp')
        
        % Remove Every Nth Trace for Data Reduction
        rmNtrc = 2;
        
        parfor (jj =  1:nChan, nWorkers)
%         for jj = 1:nChan

            % Extract Full-fold Traces & Sort Antenna Positions
            gatherLength = size(Radar{jj,ii},2); % Length of Each Channel
            
            % Flag Un-binned Traces
            if jj <= traceMod
                xTrc = 1;
            else
                xTrc = 0;
            end
                
                % Remove un-Binned Common Offset Traces
                Radar{jj,ii} = Radar{jj,ii}(:,1:gatherLength - xTrc);
                
                % Remove un-Binned Trace Indicies
                traceIx{jj,ii} = traceIx{jj,ii}(:,1:gatherLength - xTrc);
                
                % Reduce Data Volume
                if isReduceData
                    Radar{jj,ii} = Radar{jj,ii}(:,1:rmNtrc:end);
                    traceIx{jj,ii} = traceIx{jj,ii}(1:rmNtrc:end);
                end
                
                % Create Temporary Two-Way Time
                tmpTime = [0:dt:(size(Radar{jj,ii},1)-(1+padding))*dt];
                
                % Apply Systematic Channel Shifts
                if chanShift(jj) < 0
                    tmp = padarray(Radar{jj,ii},[abs(chanShift(jj)),0],...
                        mean(Radar{jj,ii}(1,:)),'pre');
                    Radar{jj, ii} = tmp(1:length(tmpTime),:);

                elseif chanShift(jj) > 0
                    tmp = padarray(Radar{jj,ii},[abs(chanShift(jj)),0],...
                        mean(Radar{jj,ii}(1,:)),'post');
                    Radar{jj,ii} = tmp(abs(chanShift(jj))+1:end-padding,:);
                   
                else
                    Radar{jj,ii} = Radar{jj,ii};
                  
                end
                
                % Process Common Offset Channels
%                 disp(' ')           
%                 fprintf(['Begin Signal Processing in Common-Offset Domain ',...
%                     filename, ' CHAN0', num2str(jj),'\n'])

                Radar{jj,ii} = processCommonOffset(Radar{jj,ii}, f0, dt );

        end
        if isReduceData
            tmpIx = sort(cat(2,traceIx{:,ii}));
            trhd{ii} = trhd{ii}(:,tmpIx);
        end
        % Store Travel-Time Axis
        TimeAxis{ii} = [0:dt:(dt.*(size(Radar{1,ii},1)-1))]';
        
        fprintf('Signal Processing Done \n')
        toc
        display(' ')
    end
    clear('Rad','tmpIx')