%% Pick Travel Time Horizons
% The Multi-channel Data Gathers are Gained for Visualization and Run
% through polarPicker.m the GUI predictive picking algorithm.
% The user must first pick the Air Wave, This is a Direct Wave! TravelTime
% Horizons must also be selected in order of increasing travel-time.
% Any number of Subsequent Direct or Reflected Wave may be chosen for this
% Analysis. T However the MxHL Product will utilize only the first Surface
% Wave and First Reflection Horizon for the Density Modeling.

% AGC Gain for PolarPicker 
    for ii = 1:nFiles
        parfor (jj = chan, nWorkers)
            directRadar{jj,ii} = AGCgain(directRadar{jj,ii},50,2);
        end
    end
    
    % Pick Direct Wave Arrival
    display('Pick Direct Wave')
    display(' ')
    display('Initial Horizon must be the Air Wave!')
    display(' ')

    [~, DirectFBpick] = polarPickerT8(directRadar);
        
    reflectRadar = Radar;
    
        % AGC Gain for PolarPicker 
    for ii = 1:nFiles
        parfor (jj = chan, nWorkers)
            reflectRadar{jj,ii} = AGCgain(reflectRadar{jj,ii},350,2);
        end
    end
    
    % Plot Common-offset gathers
    isPlotOffsetGathers = 1;
    if isPlotOffsetGathers
        for ii = 1:nFiles
            for jj = chan
                figure(chan(jj)+(ii-1).*length(chan)+100);...
                    imagesc(reflectRadar{jj,ii});colormap(LateNite);
%                 figure(chan(jj)+(ii-1).*length(chan)+200);...
%                     imagesc(directRadar{jj,ii});colormap(LateNite);
            end
        end
        
        % Pause to Examine Common-offset Gathers
        display('Review the Common-Offet Gathers')
        display('Press F5 to Continue')
        display(' ')
        keyboard
    end
    
    % Pick Primary Reflection Arrival
    display('Pick Primary Reflection')
    display(' ')
    display('Reflection Horizons Must be Picked Sequenctially in Time!')
    display(' ')
    [~, ReflectionFBpick] = polarPickerT8(reflectRadar);
    
    clear reflectRadar
    
    % Determine Number of Direct Wave Horizons
    nDirectHorizon = size(DirectFBpick,2);
    % Determine Number of Reflection Horizons
    nReflectionHorizon = size(ReflectionFBpick,2);
        
    % Convert Samples to ns
    for ii = 1:nFiles
        time{ii} = [0:dt:(size(Radar{ii,1},1)-(1+padding(ii)))*dt];
        for jj = chan
            for dh = 1:nDirectHorizon
                DirectFBpick{jj,dh,ii} = DirectFBpick{jj,dh,ii}.*dt;
            end
            for rh = 1:nReflectionHorizon
                ReflectionFBpick{jj,rh,ii} = ReflectionFBpick{jj,rh,ii}.*dt;
            end
        end
    end