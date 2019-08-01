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
            ageRadar{ii} = AGCgain(ageRadar{ii},350,2);
    end
    
    % Plot Common-offset gathers
    isPlotAgeGather = 0;
    if isPlotAgeGather
        for ii = 1:nFiles
            for jj = 1:length(RadarStack)
                figure();
                    imagesc(ageRadar{jj,ii});colormap(cmapAdapt(ageRadar{jj,ii},colorbrew));
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
    display('Isochronous Horizons Must be Picked Sequentially in Time!')
    display(' ')
    [~, isochronePick] = polarPickerT8(ageRadar);
    
    clear ageRadar isPlotAgeGather