%% Pick Depth-Horizons
% After conversion to depth, some residual between the MxHL
% Age-Depth model and the radar isochrones exists. Use this program with
% polarPicker.m to identify isochrones and calculate the modeled misfit.
% This update is applied to the isochronogram after fx-deconvolution.

% Allocation
    depthRadar = RadarDepth;
    
        % AGC Gain for PolarPicker 
    for ii = 1:nFiles
            depthRadar{ii} = AGCgain(depthRadar{ii},350,2);
    end
    
    % Plot Common-offset gathers
    isPlotDepthImage = 0;
    if isPlotDepthImage
        for ii = 1:nFiles
            for jj = 1:length(RadarStack)
                figure();
                    imagesc(depthRadar{jj,ii});colormap(cmapAdapt(depthRadar{jj,ii},colorbrew));
            end
        end
        
        % Pause to Examine Deposition Depth Image
        display('Review the Depth Image')
        display('Press F5 to Continue')
        display(' ')
        keyboard
    end
    
    % Pick Isochrones
    display('Pick Isochrones')
    display(' ')
    display('Isochronous Horizons Must be Picked Sequentially in Depth!')
    display(' ')
    [~, depthPick] = polarPickerT8(depthRadar);
    
    clear('depthRadar','isPlotDepthImage')