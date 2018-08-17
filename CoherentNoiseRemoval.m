    % Moving Median Subtraction Filter
    % This Code Applies a tapered window the Radar Data
    % The top of the time-section is unfiltered to maintain the airwave and
    % surface wave arrivals.
    
    medR = 1000; % Rank of Median Subtraction Window: W = 2R+1
    % Subtract Median at Later Times
    if f0 == 500
        MedT = xcorrWindow(2,:)+350;
        Taper = 0.25;
        if isTrimTWT
            for ii = 1:nFiles
                for jj = 1:nChan
                    if MedT(jj) > size(Radar{jj,ii},1)
                        MedT(jj) = size(Radar{jj,ii},1);
                    end
                end
            end
        end
    end
    
    if f0 == 1000
        MedT = ones(1,nChan);
        Taper = 0.25;
    end
    
    for ii = 1:nFiles
        parfor (jj = chan, nWorkers)
%         for jj = chan;
            % N samples of topTaper
            topTaper = round((MedT(jj)./2).*Taper);
            % Allocate Unfiltered Signal Matrix
            topRadar = zeros(size(Radar{jj,ii}));
            % Allocate Median Removal Matrix
            medRadar = zeros(size(Radar{jj,ii}));            
            % Solve Percentage for Median Window Tukey Taper
            medTaper = 2.*(topTaper./(size(medRadar,1)-MedT(jj)+topTaper));          
            
             topRadar(1:MedT(jj),:) = Radar{jj,ii}(1:MedT(jj),:);
             topRadar(1:MedT(jj),:) = topRadar(1:MedT(jj),:).*(tukeywin(MedT(jj),Taper)*ones(1,size(topRadar,2)));
                          
             medRadar(MedT(jj)-(topTaper-1):end,:) = movingMedianSubtraction...
                 (Radar{jj,ii}(MedT(jj)-(topTaper-1):end,:),medR);
             
             tmedL = size(medRadar,1)-MedT(jj)+topTaper;
             medRadar(MedT(jj)-(topTaper-1):end,:).*(tukeywin(tmedL,medTaper)*ones(1,size(medRadar,2)));
             Radar{jj,ii} = topRadar+medRadar;
        end
    end