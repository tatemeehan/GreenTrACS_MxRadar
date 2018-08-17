%% Read Sensors and Software Data
% This Script Reads the Binary Multi-channel GPR Data and Performs the
% pre-processing and signal processing routines. 
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
    Year = zeros(nFiles,1);
    for ii = 1 : nFiles
        tic
        %------------------------------------------------------------------
        % Multiplexed Channel Record
        if lineNo(ii)<10
            filename=['LINE0' num2str(lineNo(ii))];
        else
            filename=['LINE' num2str(lineNo(ii))];
        end
        
        % Read Data
        filepath = fullfile(dataDir,filename);
        [Rad{ii},hdr1,trhd{ii},dt,f0,~,dx,date] = readSensorsSoftwareData( filepath );
        Year(ii) = str2num(date(1:4));
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
            
            trhd{ii}(:,dupIx) = [];         % Remove Static Trace Headers from Multiplexed Record
            trhd{ii}(1,:) = 1:length(trhd{ii}); % Configure Trace Indicies
            trhd{ii}(2,:) = [0:dx:(length(trhd{ii})-1)*dx]; % Configure Distance
            Rad{ii}(:,dupIx) = [];      % Remove Static Traces from Multiplexed Data
            xArray = trhd{ii}(2,:);         % Define Configured Distance xArray
            
            fprintf('Static Trace Removal Done \n')
            toc
            display(' ')
            
        end
        
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
        
        % Pad Data with Instrument Zero
        if strcmp(date(1:4),'2016')
            padding = 50;
        end
        if strcmp(date(1:4),'2017')
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
            Array = cell(nChan,nFiles);
        end
        
        parfor (jj =  1:nChan, nWorkers)
%         for jj = chan
            % DeMux Sequential Data
            [Radar{jj,ii}, traceIx{jj,ii}, Array{jj,ii}] = DeMux(Rad{ii}, trhd{ii}, liveChan(jj));
            
            % Extract Full-fold Traces & Sort Antenna Positions
            gatherLength = length(Array{jj,ii}(1,:)); % Length of Each Chan
            traceMod = mod(length(xArray),nChan); % Un-folded Channel Index
            
            % Flag Un-binned Traces
            if jj <= traceMod
                xTrc = 1;
            else
                xTrc = 0;
            end
                
                % Remove un-Binned Common Offset Traces
                Radar{jj,ii} = Radar{jj,ii}(:,1:gatherLength - xTrc);
                
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
                disp(' ')           
                fprintf(['Begin Signal Processing in Common-Offset Domain ',...
                    filename, ' CHAN0', num2str(jj),'\n'])

                Radar{jj,ii} = processCommonOffset(Radar{jj,ii}, f0, dt );               

        end
        fprintf('Signal Processing Done \n')
        toc
        display(' ')
    end