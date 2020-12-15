function [trhd] = deadReckon( trhd, chan, delta )
	
	% Extract Multiplexed Traces
	trcIx = find(trhd(23,:)==chan);
	
	% Initialize Delta if Not Supplied
	try delta;
	catch
	    % Init 0 perturbations
	    delta = [0,0,0];
	end
	
	   % Dead Reckoning
	   R = 101; % Smoothing Window for Elevation Perturbations
	   theta = (mod(trhd(19,trcIx),90));
	           for kk = 1:length(trcIx)
	               % Convert Heading Relative to Left Axis of Quadrant
	               tmp = mod(trhd(19,trcIx(kk)),360);
	               if tmp >= 0 && tmp < 90 %Q1
	                   H = sqrt(delta(1).^2+(delta(2)-trhd(21,trcIx(kk))).^2);
	                   dl = delta(1).*sind(theta(kk));
	                   dy = abs((delta(2)-trhd(21,trcIx(kk))).*cosd(theta(kk))-dl);
	                   dx = sqrt(H.^2-dy.^2);
	                   dz = movmean([delta(2)-trhd(21,trcIx(kk))].*sind(trhd(17,trcIx(kk))),R)+delta(3);
	                   % Overwrite Antenna Midpoint Positions in Trace Header
%                        trhd(13:15,trcIx(kk)) = trhd(13:15,trcIx(kk)) + [-dx;-dy;-dz];
	                   trhd(13:15,trcIx(kk)) = trhd(13:15,trcIx(kk)) + [-dx;-dy;dz];
	               elseif tmp >= 90 && tmp < 180%Q4
	                   H = sqrt(delta(1).^2+(delta(2)-trhd(21,trcIx(kk))).^2);
	                   dl = delta(1).*sind(theta(kk));
	                   dx = (delta(2)-trhd(21,trcIx(kk))).*cosd(theta(kk))-dl;
	                   dy = sqrt(H.^2-dx.^2);
	                   dz = movmean([delta(2)-trhd(21,trcIx(kk))].*sind(trhd(17,trcIx(kk))),R)+delta(3);
	                   % Overwrite Antenna Midpoint Positions in Trace Header
%                        trhd(13:15,trcIx(kk)) = trhd(13:15,trcIx(kk)) + [-dx;dy;-dz];
	                   trhd(13:15,trcIx(kk)) = trhd(13:15,trcIx(kk)) + [dx;dy;dz];
	               elseif tmp >= 180 && tmp < 270%Q3
	                   H = sqrt(delta(1).^2+(delta(2)-trhd(21,trcIx(kk))).^2);
	                   dlx = delta(1).*cosd(theta(kk));
	                   dly = delta(1).*sind(theta(kk));
	                   dx = abs((delta(2)-trhd(21,trcIx(kk))).*sind(theta(kk))+dlx);
	                   dy = abs((delta(2)-trhd(21,trcIx(kk))).*cosd(theta(kk))-dly);
	                   dz = movmean([delta(2)-trhd(21,trcIx(kk))].*sind(trhd(17,trcIx(kk))),R)+delta(3);
	                   % Overwrite Antenna Midpoint Positions in Trace Header
                       trhd(13:15,trcIx(kk)) = trhd(13:15,trcIx(kk)) + [dx;dy;dz];
% 	                   trhd(13:15,trcIx(kk)) = trhd(13:15,trcIx(kk)) + [dx;dy;-dz];
	               else%Q2
	                   H = sqrt(delta(1).^2+(delta(2)-trhd(21,trcIx(kk))).^2);
	                   dl = delta(1).*sind(theta(kk));
	                   dx = abs((delta(2)-trhd(21,trcIx(kk))).*cosd(theta(kk))-dl);
	                   dy = sqrt(H.^2-dx.^2);
	                   dz = movmean([delta(2)-trhd(21,trcIx(kk))].*sind(trhd(17,trcIx(kk))),R)+delta(3);
	                   % Overwrite Antenna Midpoint Positions in Trace Header
%                        trhd(13:15,trcIx(kk)) = trhd(13:15,trcIx(kk)) + [dx;-dy;-dz];
	                   trhd(13:15,trcIx(kk)) = trhd(13:15,trcIx(kk)) + [dx;-dy;dz];
	               end
	           end
	
	            % Distance Array Corrected to MidPoint Centers
	            trhd(16,trcIx) = trhd(16,trcIx) + trhd(21,trcIx);
	            % Unwrap Heading to 0 - 360
	            trhd(19,trcIx) = mod(trhd(19,trcIx),360);
	            % Ensure that no points equal 360 exactly
	            ix360 = trhd(19,trcIx)==360;
	            trhd(19,trcIx(ix360))=0;
	end

