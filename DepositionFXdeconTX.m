    if isDepthFXdecon
    % Spatial Predictive Filtering for Random Noise Suression
    Window = 2/dt; % Number of Temporal Samples in Fx-Decon Window Per Computation
    Overlap = round(Window/2); % Number of Overlapping Computations
    mu = 1e5; % Weighting
    nfR = 50; % Length of Deconvolution Filter
    fxProRadar = cell(nFiles,1); % Allocation
    
    for ii = 1:nFiles
%         parfor (jj = chan, nWorkers)
%         for jj = chan
            tmp = 1:length(RadarDeposition{ii}(:,1));
            filterLap = Overlap:Overlap:length(tmp);
            % Design Data Window
            for kk = filterLap
                if kk == filterLap(1)
                    run = 1;
                    traceIx = 1:kk+(Overlap);
                    tmpData = RadarDeposition{ii}(traceIx,:);
                elseif kk == filterLap(end)
                    run = length(filterLap);
                    traceIx = kk-(Overlap):length(tmp);
                    tmpData = RadarDeposition{ii}(traceIx,:);
                else
                    run = find(kk==filterLap);
                    traceIx = kk-(Overlap):kk+(Overlap);
                    tmpData = RadarDeposition{ii}(traceIx,:);
                end
                % Perform Fx-Decon within Windowed Data
                fxProRadar{ii}{run,1} = fx_deconT8(tmpData,dt,nfR,mu,f0);
                fxProRadar{ii}{run,2} = traceIx;
%             end
        end
        
        fxRadar = cell(nFiles,1);   % Allocate Temporary Arrays
        fxCatRadar = cell(nFiles,1);% Temporary Array
        
%         parfor (jj = chan, nWorkers)
            % Concatenate Deconvolved Matricies
            fxCatRadar{ii}{1} = cat(1,fxProRadar{ii}{:,1});
            % Sort Trace Indicies
            [fxCatRadar{ii}{2},fxCatRadar{ii}{3}] = sort([fxProRadar{ii}{:,2}]);
            fxRadar{ii} = fxCatRadar{ii}{1}(fxCatRadar{ii}{3},:);

            % Fold Replicated Traces by Mean Stacking
            for kk = 1:length(RadarDeposition{ii}(:,1));
                RadarDeposition{ii}(kk,:) = mean(fxRadar{ii}(find(fxCatRadar{ii}{2}==kk),:),1);
            end
%         end
    end
    clear('fxCatRadar','fxProRadar','fxRadar','stackIx')
    end
    
    if isSpatialFXdecon
        
        Window = 50; % Number of Traces in Fx-Decon Window Per Computation
        Wn = 2;      % Number of Overlapping Computations
        Overlap = round(Window/Wn);
        mu = 1e3;   % Weighting
        nfR = 25; % Length of Deconvolution Filter
        fxProRadar = cell(nFiles,1); % Allocation
        
        for ii = 1:nFiles
%             for jj = 1
%             parfor (jj = chan, nWorkers)
                tmpfx = 1:length(RadarDeposition{jj,ii});
                filterLap = Overlap:Overlap:length(tmpfx);
                % Design Data Window
                for kk = filterLap
                    if kk == filterLap(1)
                        run = 1;
                        traceIx = 1:kk+(Overlap);
                        tmpData = RadarDeposition{jj,ii}(:,traceIx);
                    elseif kk == filterLap(end)
                        run = length(filterLap);
                        traceIx = kk-(Overlap):length(tmpfx);
                        tmpData = RadarDeposition{jj,ii}(:,traceIx);
                    else
                        run = find(kk==filterLap);
                        traceIx = kk-(Overlap):kk+(Overlap);
                        tmpData = RadarDeposition{jj,ii}(:,traceIx);
                    end
                    % Perform Fx-Decon within Windowed Data
                    fxProRadar{jj,ii}{run,1} = fx_deconT8(tmpData,dt,nfR,mu,f0);
                    fxProRadar{jj,ii}{run,2} = traceIx;
                end
%             end
            
            fxRadar = cell(nFiles,1);   % Allocate Temporary Arrays
            fxCatRadar = cell(nFiles,1);% Temporary Array
            
%             parfor (jj = chan, nWorkers)
%             for jj = 1;
                % Concatenate Deconvolved Matricies
                fxCatRadar{ii}{1} = cat(2,fxProRadar{ii}{:,1});
                % Sort Trace Indicies
                [fxCatRadar{ii}{2},fxCatRadar{ii}{3}] = sort([fxProRadar{ii}{:,2}]);
                fxRadar{ii} = fxCatRadar{ii}{1}(:,fxCatRadar{ii}{3});
                
                % Fold Replicated Traces by Mean Stacking
                for kk = 1:length(RadarDeposition{ii}(1,:));
                    RadarDeposition{ii}(:,kk) = mean(fxRadar{ii}(:,find(fxCatRadar{ii}{2}==kk)),2);
                end
%             end

        end
        clear('fxCatRadar','fxProRadar','fxRadar','out','stackIx')
    end
       % Normalize Image
        RadarDeposition{ii} = trcNormalize(RadarDeposition{ii});
