%% Pick Age-Horizons
% After conversion to depositional stratum, some residual between the MxHL
% Age-Depth model and the radar isochrones exists. Use this program with
% polarPicker.m to identify isochrones and calculate the modeled misfit.
% This update is applied to the isochronogram prior to fx-deconvolution.

% Allocation
    ageRadar = RadarDeposition;
    deltaAgeModel = cell(nFiles,1);
    depositionAgeModel = cell(nFiles,1);
    
        % AGC Gain for PolarPicker 
    for ii = 1:nFiles
%         parfor (jj = chan, nWorkers)
            ageRadar{ii} = AGCgain(ageRadar{ii},350,2);
%         end
    end
    
    % Plot Common-offset gathers
    isPlotOffsetGathers = 0;
    if isPlotOffsetGathers
        for ii = 1:nFiles
            for jj = chan
                figure();
                    imagesc(ageRadar{jj,ii});colormap(colorbrew);
%                 figure(chan(jj)+(ii-1).*length(chan)+200);...
%                     imagesc(directRadar{jj,ii});colormap(LateNite);
            end
        end
        
        % Pause to Examine Deposition Age Gather
        display('Review the Deposition Age Gather')
        display('Press F5 to Continue')
        display(' ')
        keyboard
    end
    
    % Pick Isochrones
    display('Pick Isochrones')
    display(' ')
    display('Isochronous Horizons Must be Picked Sequenctially in Time!')
    display(' ')
    [~, isochronePick] = polarPickerT8(ageRadar);
    
    clear ageRadar
    
    % Determine Number of Age Horizons
    nAgeHorizon = size(isochronePick,2);
    
    % Calculate Residual
    ageResidual = isochronePick;
    datumAge = zeros(nAgeHorizon,nFiles);
    warning('off','MATLAB:mir_warning_changing_try_catch');
    try iceCoreIx;
    catch
        iceCoreIx = 1:100;
    end
        
    for ii = 1:nFiles
        dAge = mean(diff(DepositionAxis{ii}));
        for jj = 1;% chan (Loop over chan for Mx processing)
            for ah = 1:nAgeHorizon
                % Convert Samples to Age
                isochronePick{jj,ah,ii} = isochronePick{jj,ah,ii}.*dAge;
                % Determine Measurement Datum
                [~,datumIx] = min(abs(mean(isochronePick{jj,ah,ii}(iceCoreIx))-DepositionAxis{ii}));
                datumAge(ah,ii) = DepositionAxis{ii}(datumIx);
                % Calculate Misfit
                ageResidual{jj,ah,ii} = isochronePick{jj,ah,ii} - datumAge(ah,ii);
            end
        end
    end