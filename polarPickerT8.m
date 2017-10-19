function [pickX,pickT] = polarPickerT8( Radar, snap, nHorizon, isZoom )
% polarPicker is a predicitive picking algorithm which applies a polar
% transform to a window of the radargram and seeks the direction theta for
% which the specified amplitude gradient is least. The window advances
% 'step' traces and recomputes theta for the next pick prediction. The seed
% is manually selected by the geophysicsist and a plot displays the picks 
% per iteration. Raw Picks are smoothly interpolated for all X indicies in
% Radar Gram. The Look-ahead distance and foreward step distance were found
% to be optimum at 3 and 25 and have been hardcoded. 
% polarPicker.m will snap picks to a peak, trough, or zero-crossing and
% handles multiple horizons for selection. polarPicker.m incorporates a GUI
% which allows the geophysicsist to 'Re-Seed' the pick if the intial 
% conditioning goes awry. Additionally, if snap and nHorizon are 
% unspecified the user will be propmted via the GUI to supply this
% information.
%
%   The idea for this picking routine is founded upon the mountain range
%   hiker metaphor by Alex Miller. A hiker will follow the path requiring
%   the least work to cross the Radar Wiggle Mountain Range.
%
%   Inputs: Radar - Cell-Array of Processed Radar Channels and Files 
%           Optional:
%               snap  - string indicating pick snap: 'peak','trough','zero'
%               nHorizon - Interger Indicated the Number of Horzons to Pick
%               isZoom - Toggle 0/1 for Zoom Trace Display
%
%   Output: pickX - Cell-Array of interpolated X indicies
%           pickT - Cell-Array of interpolated T indicies
%
%   Written by Tate Meehan, Boise State University, GreenTracs 2017

% Surpress Rounding Warning in Integer Indexing Operations
ID = 'MATLAB:colon:nonIntegerIndex';
warning('off', ID);

% Load Custom colormap
cmap = load('LateNite.mat');
cmap = struct2cell(cmap);
cmap = cell2mat(cmap);
% Establish Mask Color for Picks
pColorIx = round(quantile(1:length(cmap),.5));
tColorIx = round(quantile(1:length(cmap),.5));
peakColor = cmap(pColorIx,:);
troughColor = cmap(tColorIx,:);
zeroColor = 1-mean([peakColor;troughColor],1);

% Determine Cell Array Dimensions for looping
[nChan,nFiles] = size(Radar);

% Get User Input for nHorizon: Default is 1
if nargin <= 2
    horizonQuest = inputdlg('  Number of Horizons to Pick  :',...
        'Polar Picker', 1,{'1'});
    if isempty(horizonQuest);
        % Horizons is 0
        nHorizon = 0;
    else
        nHorizon = str2num(horizonQuest{:});
    end
end

% Set Zoom Toggle: Default is On
if nargin <= 3
    isZoom = 1;
end

% Allocate Cell Arrays for Pick Export
pickX = cell(nChan,nHorizon,nFiles);
pickT = cell(nChan,nHorizon,nFiles);
pickTT = cell(nChan,nHorizon,nFiles);
pickTP = cell(nChan,nHorizon,nFiles);

% Create Re-Seed Callback
global isPicking
isPicking = 1;

if nHorizon > 0
for ii = 1:nFiles
    jj = 1;
    while jj <= nChan
        kk = 1;
        % Open RadarGram for Picks
        figstr = ['Polar Picker: Channel: ',num2str(jj,'%02d'),' File: ' num2str(ii,'%02d')];
        h = figure('Name',figstr,'NumberTitle','off');
%         set(h,'units','normalized','outerposition',[.25 .25 .5 .5 ])
        set(h,'units','normalized','outerposition',[.25 .25 .25 1 ])
%         set(h,'units','normalized','outerposition',[0 0 1 1 ])
        btn = uicontrol('Style', 'pushbutton', 'String', 'Re-Seed',...
           'FontWeight','bold','FontSize',12,'Position', [125 10 100 50],...
            'Callback',{@reSeed});
        imagesc(Radar{jj,ii});colormap(cmap);hold on;
        
        while kk <= nHorizon
            zoom out;   % Refresh orignal trace display
        % Construct a questdlg to Choose Pick Snap
        if nargin < 2;
            prompt = ['*',blanks(2),'Which wavelet phase will you pick?',blanks(2),'*'];
            titleTxt = ['File: ' num2str(ii,'%02d'),' Channel: ',...
                num2str(jj,'%02d'),' Horizon: ',num2str(kk,'%02d')];
            snapQuest = questdlg(prompt, titleTxt,'Peak','Trough','Zero','Peak');
            if isempty(snapQuest)
                % Default is Peak
                snapQuest = 'Peak';
            end
            if strcmp(snapQuest,'Peak')
                snap = 'peak';
            elseif strcmp(snapQuest,'Trough')
                snap = 'trough';
            elseif strcmp(snapQuest,'Zero')
                snap ='zero';
            end
        end
        % Zoom Trace Dispay
        if isZoom
            zoom on;
            display(['Zoom Trace Display: Press Any Key to Continue'])
            display(' ')
            pause()
            zoom off;
        end
        
        [x,t] = ginput(1); % Plant the Seed
        x = floor(x); t = floor(t); % Round Index to Integer Value
        if x == 0;
            x = 1;  % Check Array Index
        end
        
        look = 3;   % Window Dimension (2.*look+1,look)
        fore = 10;  % Step Distance for Next Prediction
        pace = 1;   % Counter
        
        % Prevent Bad Seed Pick and Subsequent Errors
        while x < 1 || x > size(Radar{jj,ii},2) || t < 1 || ...
                t > size(Radar{jj,ii},1) || (size(Radar{jj,ii},2))-x <= look-fore
            display('Pick Index Out of Range: Re-Pick the Seed')
            [x,t] = ginput(1); % Plant the Seed
            x = floor(x); t = floor(t); % Round Index to Integer Value
            if x == 0;
                x = 1;  % Check Array Index
            end
        end
        
        pickT{jj,kk,ii}(pace) = t;
        pickX{jj,kk,ii}(pace) = x;
        
        % Plot Seed Pick on Current Figure
        if strcmp(snap,'peak')
            figure(gcf);plot(x,t,'.','color',peakColor);
        end
        if strcmp(snap,'trough')
            figure(gcf);plot(x,t,'.','color',troughColor);
        end
        if strcmp(snap,'zero')
            tp = t;
            tt = t;
            pickTP{jj,kk,ii}(pace) = tp;
            pickTT{jj,kk,ii}(pace) = tt;
            figure(gcf);plot(x,tt,'.','color',troughColor);
            figure(gcf);plot(x,tp,'.','color',peakColor);
        end
        
        % Polar Coordinate Transform
        % Matrix of T Coodinate Indicies Forward Centered on Pick
        T = repmat(-look:look,look,1)';
        % Matrix of X Coodinate Indicies Forward Centered on Pick
        X = repmat(0:look-1,2.*look+1,1);
        % Matrix of Radial Distances Forward Centered on Pick
        D = sqrt(T.^2+X.^2);
        % Matrix of Degree Theta Forward Centered on Pick
        O = asind(X./D);
        % Find NaN at Origin
        nanIx = find(isnan(O));
        [NaNjj,NaNii] = ind2sub(size(O),nanIx);
        % Replace NaN with 90 Degrees
        O(NaNjj,NaNii) = 90;
        % Keep Theta Sign Information in Quadrant 2
        O(1:look,:) = 180-O(1:look,:);
        
        while x<size(Radar{jj,ii},2)-look
            % Check CallBack Command for RePick
            if ~isPicking
                [x,t] = ginput(1); % Plant the Seed
                x = floor(x); t = floor(t); % Round Index to Integer Value
                if x == 0;
                    x = 1;  % Check Array Index
                end
                
                % Prevent Bad Seed Pick and Subsequent Errors
                while x < 1 || x > size(Radar{jj,ii},2) || t < 1 || ...
                        t > size(Radar{jj,ii},1) || (size(Radar{jj,ii},2))-x <= look-fore
                    display('Pick Index Out of Range: Re-Pick the Seed')
                    display(' ')
                    [x,t] = ginput(1); % Plant the Seed
                    x = floor(x); t = floor(t); % Round Index to Integer Value
                    if x == 0;
                        x = 1;  % Check Array Index
                    end
                end
                               
                % Replace Pick Index with RePickIx
                tmp = pace; % Bookmark Old Index
                [~, pace] = min(abs(pickX{jj,kk,ii}-x));
                pickT{jj,kk,ii}(pace) = t;
                pickX{jj,kk,ii}(pace) = x;
                % Clear Bad Picks from Figure
                items = get(gca,'Children');
                nBadPicks = tmp - pace + 1;
                % If User Skips Ahead Remove No Picks
                if nBadPicks <= 0
                    nBadPicks = [];
                end
                if strcmp(snap,'zero')
                    % Remove Peak and Trough Picks if Zero-Corssing
                    nBadPicks = 2.*nBadPicks;
                end
                % Delete Graphics
                delete(items(1:nBadPicks))
                % Clear Graphics Memory
                items(1:nBadPicks) = [];
                
                % Plot ReSeed Pick on Current Figure
                if strcmp(snap,'peak')
                    figure(gcf);plot(x,t,'.','color',peakColor);
                end
                if strcmp(snap,'trough')
                    figure(gcf);plot(x,t,'.','color',troughColor);
                end
                if strcmp(snap,'zero')
                    tp = t;
                    tt = t;
                    pickTP{jj,kk,ii}(pace) = tp;
                    pickTT{jj,kk,ii}(pace) = tt;
                    figure(gcf);plot(x,tt,'.','color',troughColor);
                    figure(gcf);plot(x,tp,'.','color',peakColor);
                end
                
                % Switch Picking Back on
                isPicking = 1;
            end
            
            if strcmp(snap,'zero')
                % Condition Evaluates Peak and Trough Prediction if 'zero'
                % Run Peak Picker
                % Grab Appropriate Data Window
                if tp-look<1
                    trailP = Radar{jj,ii}(1:tp+look,x:x+look-1);
                    Otheta = O(((2*look+1)-size(trailP,1))+1:(2*look+1),:);
                    Dradius = D(((2*look+1)-size(trailP,1))+1:(2*look+1),:);
                elseif tp+look>size(Radar{jj,ii},1)
                    trailP = Radar{jj,ii}(tp-look:size(Radar{jj,ii},1),x:x+look-1);
                    Otheta = O(1:size(trailP,1),:);
                    Dradius = D(1:size(trailP,1),:);
                else
                    trailP = Radar{jj,ii}(tp-look:tp+look,x:x+look-1);
                    Otheta = O;
                    Dradius = D;
                end
                % Normalize Windowed Data
                if abs(min(trailP(:))) > max(trailP(:))
                    maxValue = abs(min(trailP(:)));
                    minValue = min(trailP(:));
                else
                    maxValue = max(trailP(:));
                    minValue = -max(trailP(:));
                end
                trailP = 2 .* trailP ./ (maxValue - minValue);
                                                
                % Concatenate Array of Peak Theta, Distance, Amplitude
                MapP = [Otheta(:),Dradius(:),trailP(:)];
                
                % Run Trough Picker
                if tt-look<1
                    trailT = Radar{jj,ii}(1:tt+look,x:x+look-1);
                    Otheta = O(((2*look+1)-size(trailT,1))+1:(2*look+1),:);
                    Dradius = D(((2*look+1)-size(trailT,1))+1:(2*look+1),:);                    
                elseif tt+look>size(Radar{jj,ii},1)
                    trailT = Radar{jj,ii}(tt-look:size(Radar{jj,ii},1),x:x+look-1);
                    Otheta = O(1:size(trailT,1),:);
                    Dradius = D(1:size(trailT,1),:);                    
                else
                    trailT = Radar{jj,ii}(tt-look:tt+look,x:x+look-1);
                    Otheta = O;
                    Dradius = D;                    
                end
                % Normalize Windowed Data
                if abs(min(trailT(:))) > max(trailT(:))
                    maxValue = abs(min(trailT(:)));
                    minValue = min(trailT(:));
                else
                    maxValue = max(trailT(:));
                    minValue = -max(trailT(:));
                end
                trailT = 2 .* trailT ./ (maxValue - minValue);
                
                % Concatenate Array of Peak Theta, Distance, Amplitude
                MapT = [Otheta(:),Dradius(:),trailT(:)];
                
                clear ('trailT','trailP')
            else
                % Run Picker for 'peak' or 'trough' snap
                if t-look<1
                    trail = Radar{jj,ii}(1:t+look,x:x+look-1);
                    Otheta = O(((2*look+1)-size(trail,1))+1:(2*look+1),:);
                    Dradius = D(((2*look+1)-size(trail,1))+1:(2*look+1),:);                        
                elseif t+look>size(Radar{jj,ii},1)
                    trail = Radar{jj,ii}(t-look:size(Radar{jj,ii},1),x:x+look-1);
                    Otheta = O(1:size(trail,1),:);
                    Dradius = D(1:size(trail,1),:);                       
                else
                    trail = Radar{jj,ii}(t-look:t+look,x:x+look-1);
                    Otheta = O;
                    Dradius = D;                        
                end
                
                % Normalize Windowed Data
                if abs(min(trail(:))) > max(trail(:))
                    maxValue = abs(min(trail(:)));
                    minValue = min(trail(:));
                else
                    maxValue = max(trail(:));
                    minValue = -max(trail(:));
                end
                trail = 2 .* trail ./ (maxValue - minValue);
                               
                % Concatenate Array of Peak Theta, Distance, Amplitude
                Map = [Otheta(:),Dradius(:),trail(:)];
                
                clear ('trail','map')
            end
            
            % Using the Hiker Analogy - The polar space transform
            % optimizes the location where the mountain gradient is
            % steepest. A peak is represented by a maximum where 
            % work opposing gravity will be the greatest for the hiker.
            % A trough in the RadarGram represents the minimum value
            % in the polar space cost function. The Zero contour is taken
            % to be the average of the maximum and minimum cost diections.
                                  
            if strcmp(snap,'peak')
                % Bin Theta Vaules
                [~,theta,thetaBin] = histcounts(Map(:,1),unique(Map(:,1)));
                % Grab Peak Theta Index
                [~,thetaIx] = max(accumarray(thetaBin,Map(:,3),[],@mean));
                % Grab the Map Index
                MapIx = find(thetaBin == thetaIx);
                % Grab Corresponding Distance
                distance = mean(Map(MapIx,2));
                % Compute Corresponding Time Coordinate
                predictT = distance.*cosd(theta(thetaIx));
                % Advance t with prediciton
                t = t + predictT;
                % Advance x with hard step
                x = x + fore;
                % Advance counter
                pace = pace+1;
                % Store Picks
                pickX{jj,kk,ii}(pace) = x;
                pickT{jj,kk,ii}(pace) = t;
                % Plot Current Pick
                pause on; pause(0.001);
                figure(gcf);plot(x,t,'.','color',peakColor)
                
            elseif strcmp(snap,'trough')
                % Bin Theta Vaules
                [~,theta,thetaBin] = histcounts(Map(:,1),unique(Map(:,1)));
                % Grab Trough Theta Index
                [~,thetaIx] = min(accumarray(thetaBin,Map(:,3),[],@mean));
                % Grab the Map Index
                MapIx = find(thetaBin == thetaIx);
                % Grab Corresponding Distance
                distance = mean(Map(MapIx,2));
                % Compute Corresponding Time Coordinate
                predictT = distance.*cosd(theta(thetaIx));
                % Advance t with prediciton
                t = t + predictT;
                % Advance x with hard step
                x = x + fore;
                % Advance counter
                pace = pace+1;
                % Store Pick
                pickX{jj,kk,ii}(pace) = x;
                pickT{jj,kk,ii}(pace) = t;
                % Plot Current Pick
                pause on; pause(0.001);
                figure(gcf);plot(x,t,'.','color',troughColor)
                clear ('Map')
                
            elseif strcmp(snap,'zero')
                % Bin Theta Vaules
%                 [~,theta,thetaBin] = histcounts(MapP(:,1),2.*look+1);
                [~,theta,thetaBin] = histcounts(MapP(:,1),unique(MapP(:,1)));
                % Grab Peak Theta Index
                [~,thetaIx] = max(accumarray(thetaBin,MapP(:,3),[],@mean));
                % Grab the Map Index
                MapIx = find(thetaBin == thetaIx);
                % Grab Corresponding Distance
                distanceP = mean(MapP(MapIx,2));
                % Compute Corresponding Time Coordinate
                peakT = (distanceP.*cosd(theta(thetaIx)));
                % Advance t with prediciton
                tp = tp + peakT;
                
                % Bin Theta Vaules
%                 [~,theta,thetaBin] = histcounts(MapT(:,1),2.*look+1);
                [~,theta,thetaBin] = histcounts(MapT(:,1),unique(MapT(:,1)));
                % Grab Trough Theta Index
                [~,thetaIx] = min(accumarray(thetaBin,MapT(:,3),[],@mean));
                % Grab the Map Index
                MapIx = find(thetaBin == thetaIx);
                % Grab Corresponding Distance
                distanceT = mean(MapT(MapIx,2));
                % Compute Corresponding Time Coordinate
                troughT = (distanceT.*cosd(theta(thetaIx)));
                % Advance t with prediciton
                tt = tt + troughT;
                % Advance x with hard step
                x = x + fore;
                % Advance counter
                pace = pace+1;
                % Store Picks
                pickX{jj,kk,ii}(pace) = x;
                pickTP{jj,kk,ii}(pace) = tp;
                pickTT{jj,kk,ii}(pace) = tt;
                % Plot Current Picks
                pause on; pause(0.001);
                figure(gcf);plot(x,tp,'.','color',peakColor);
                plot(x,tt,'.','color',troughColor)
                    clear ('MapT','MapP')
            end
        end
        
        % Construct a questdlg to Accept or Reject Picks
        prompt = ['*',blanks(10),'Are these picks correct?',blanks(11),'*'];
        titleTxt = ['File: ' num2str(ii,'%02d'),' Channel: ',...
            num2str(jj,'%02d'),' Horizon: ',num2str(kk,'%02d')];
        questPick = questdlg(prompt, titleTxt,'Yes','No','Yes');
        switch questPick
            case 'Yes'
            if strcmp(snap,'zero')
                % Zero-Cross is Average of Peak and Trough Picks
                pickT{jj,kk,ii} = floor(mean([pickTP{jj,kk,ii};pickTT{jj,kk,ii}],1));
            end
            
            % Clear Many Children of Picks
            Save = get(gca,'Children');
            nSavePicks = length(pickT{jj,kk,ii});
            if strcmp(snap,'zero')
                % Remove Peak and Trough Picks if Zero-Corssing
                nSavePicks = length(pickTT{jj,kk,ii})+...
                    length(pickTP{jj,kk,ii});
            end
            % Delete Graphics
            delete(Save(1:nSavePicks))
            % Clear Graphics Memory
            Save(1:nSavePicks) = [];
            
            % Smoothly Interpolate Picks Across all Traces
            pickT{jj,kk,ii} = round(nonParametricSmooth...
                (pickX{jj,kk,ii},pickT{jj,kk,ii},1:size(Radar{jj,ii},2),fore));
            % Extrapolate NaN Values on Edges
            pickT{jj,kk,ii} = inpaint_nans(pickT{jj,kk,ii},1);
            pickX{jj,kk,ii} = 1:size(Radar{jj,ii},2);
            
            % Plot Smoothed Picks
            if strcmp(snap,'peak')
                figure(gcf);
                plot(pickX{jj,kk,ii},pickT{jj,kk,ii},'.','color',peakColor)
            end
            
            if strcmp(snap,'trough')
                figure(gcf);
                plot(pickX{jj,kk,ii},pickT{jj,kk,ii},'.','color',troughColor)
            end
            
            if strcmp(snap,'zero')
                figure(gcf);
                plot(pickX{jj,kk,ii},pickT{jj,kk,ii},'.','color',zeroColor)
            end
            
            % Test Horizon Count for Iteration
            if kk == nHorizon
                % Close Picker Figure
                close(figure(h.Number))
                % Advance RadarGram
                jj = jj + 1;
            end
            % Advance Horizon
            kk = kk + 1;
            
            % Ensure isPicking
            isPicking = 1;
            
            case 'No'
            display(['RePick File: ' num2str(ii,'%02d'),' Channel: ',...
                num2str(jj,'%02d'),' Horizon: ',num2str(kk,'%02d')])
            display(' ')
            
            % Clear Bad Picks
            Reject = get(gca,'Children');
            nRejectPicks = length(pickT{jj,kk,ii});
            if strcmp(snap,'zero')
                % Remove Peak and Trough Picks if Zero-Corssing
                nRejectPicks = length(pickTT{jj,kk,ii})+...
                    length(pickTP{jj,kk,ii});
            end
            
            % Delete Graphics
            delete(Reject(1:nRejectPicks))
            % Clear Graphics Memory
            Reject(1:nRejectPicks) = [];
            % Clear Bad Picks
            pickTP{jj,kk,ii} = [];
            pickTT{jj,kk,ii} = [];
            pickT{jj,kk,ii} = [];
            pickX{jj,kk,ii} = [];
            
            % Ensure isPicking
            isPicking = 1;
            
        end
        end
    end
end
end
end
%% CallBack Function for polarPicker.m
function reSeed(src,evnt)
global isPicking
    isPicking = 0;
end
