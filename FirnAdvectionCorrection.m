%% Lateral Firn Advection Correction
% The age-depth model is regridded using new (X,Z) coordinates that are
% determined by ice sheet surface velocities that are propagated backwards
% in time. The NASA MEaSUREs Multi-year Greenland Ice Sheet Velocity Mosaic
% Version 1 data are incorporated to make the advection model.

if isLoadGPS && isMEaSUREs
    for ii = 1:nFiles
        measuresDir = '/home/tatemeehan/GreenTracs2017/MEaSUREs';
        filenameX = 'greenland_vel_mosaic250_vx_v1.tif';
        filenameY = 'greenland_vel_mosaic250_vy_v1.tif';
        surfV = surfaceVelocity(measuresDir,filenameX,filenameY,trhd{ii});
        % Compute X-Position Perturbation for Age Model
        deltaX = surfV'.*AgeModel{ii};
        % Filter Bad Perturbations due to Array Turning
        % Forward Difference
        fdiff = deltaX(end,1:end-99) - deltaX(end,100:end);
        fix = find(fdiff > quantile(fdiff,0.99) | fdiff < quantile(fdiff,0.01)) + 99;
        goodIx = find(~ismember(1:length(deltaX(end,:)),fix));
        % Filter Discontinuities by Replacement with Nearest Neighbor
        for kk = 1:length(fix)
            [~,binix] = min(abs(fix(kk)-goodIx));
            deltaX(:,fix(kk)) =  deltaX(:,binix);
        end
        
        % Correct Depth-Age Model for Lateral Firn Advection
        % Backward Propagation of X perturbations
        dxStack = double(xStack{ii}-deltaX);
        staticAgeModel = griddata(dxStack,DepthMatrix{ii},...
            updateAgeModel{ii},xStack{ii},DepthMatrix{ii},'linear');
        staticAgeModel = inpaint_nans(staticAgeModel,1);
        staticAgeModel(1,:) = 0;
        smoothStaticAgeModel = imgaussfilt(staticAgeModel,[1 25]);
        accumulationFirnUpdate = ((DepthMatrix{ii}.*AvgDensityModel{ii})./(staticAgeModel+eps));
        % Currenlty No Working Variables are Assigned
    end
end