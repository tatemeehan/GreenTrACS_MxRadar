function [ Rad, trhd, trcIx, xArray ] = DeMux( Rad, trhd, chan,isLoadGPS )
% DeMux - Demultiplexes Sensors and Software Multi-offset Radar Data.
% DeadReckon - Calculates Truer GPR Antenna Positions by High Accuracy GPS. 
  % Assumtions: The GPS Antenna is located at Tx1
  %             The GPR Array is towed straight
  %             The Shot Gather Bin Center is the Array mean
  % Errors caused by these assumtions are significant, however, inevitable.
  % The dead reckon approach is utilized to distribute the GPS position to
  % each GPR antenna.
  
% Boise State University: Tate Meehan, NASA ISGC 2019

    % DeMux Data
        trcIx = find(trhd(23,:)==chan);
        Rad=Rad(:,trcIx);
        xArray = trhd(2,trcIx);
        
%         if isLoadGPS
%             % Dead Reckoning for Antenna Midpoint Position
%             dx = trhd(21,trcIx).*sind(trhd(20,trcIx));
%             dy = trhd(21,trcIx).*cosd(trhd(20,trcIx));
%             dz = trhd(21,trcIx).*sind(trhd(17,trcIx));
%             % Overwrite Antenna Midpoint Positions in Trace Header
%             trhd(13:15,trcIx) = trhd(13:15,trcIx) + [dx;dy;dz];
%             % Distance Array Corrected to MidPoint Centers
%             trhd(16,trcIx) = trhd(16,trcIx) - trhd(21,trcIx);
%         end


end

