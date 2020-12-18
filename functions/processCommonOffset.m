function [ Rad ] = processCommonOffset( Rad, f0, dt, offset )
% processCommonOffset is a data processing Subroutine
%   
%   Written by Tate Meehan, Boise State University, GreenTrACS 2016
%
% Process Data
% Filter Parameters are Established within Function
isDisplay = 0;       % Control Text Display to Screen
isNormalize = 0;     % Control Flag to Normalize Data
isDeWoW = 0;         % Control Flag to De-Wow Data
isDeTrend = 0;       % Control Flag to De-Trend Data
isMedFilt = 1;       % Control Flag for Median Subtraction Filter
isSVDfilt = 0;       % Control Flag for SVD Component Subtraction Filter
isBandPass = 1;      % Control Flag to Band-Pass Filter Data
isTimeZero = 0;      % Control Flag for Time-Zero Correction
isExpGain = 1;       % Control Flag for Ramped Gain of Data
isACGGain = 0;       % Control Fla for AGC Gain of Data
isCleanNaN = 0;      % Control Flag for NaN Clean-up of Data
isStak = 0;          % Control Flag to Stack Data

    %----------------------------------------------------------------------      
    % Convert Units
    f0Hz = f0 * 1e6;        % [Hz]
    dtSec = dt * 1e-9;      % [s]
    [nsamp, ntrcs] = size(Rad);

    %----------------------------------------------------------------------
    % Normalize Data
    if isNormalize
        if isDisplay
        display( 'Begin Normalize')
        tic
        end
        
        Rad = normalize_dylan( Rad );
        
        if isDisplay
        display( 'Normalize Done')
        toc
        display(' ')
        end
    end    
    
    %----------------------------------------------------------------------      
    % De-WoW
    if isDeWoW
        if isDisplay
        display( 'Begin De-WOW Filter')
        tic
        end
        
        Rad = deWoW(Rad,dt,f0);
        
        if isDisplay
        display( 'De-WOW Filter Done')
        toc
        display(' ')
        end
    end
    
    %----------------------------------------------------------------------    
    % Remove Trends (De-WoW)
    if isDeTrend
        if isDisplay
        display( 'Begin De-Trend')
        tic
        end
        
        % Parameters
        isMeanFilt = 1;
        isDeStripe = 0;
        
        % Function
        if isMeanFilt
            Rad = detrend( Rad, 'constant' ); % Subtract Mean of Trace
        end
        if isDeStripe
            Rad = detrend( Rad, 'linear' ); % Subtract Linear Trend of Trace
        end
        
        if isDisplay
        display( 'De-Trend Done')
        toc
        display(' ')
        end
    end
        
    %----------------------------------------------------------------------
    % Median Subtraction Filter
    if isMedFilt
        if isDisplay
        display( 'Begin Median Subtraction Filter')
        tic
        end
        
        % Parameters
        % Rank of Median Subtraction Filter
        % Nominal Frequency Pass
        MedFiltR = 2.*(ceil(1/((f0Hz)*dtSec))-1)+1;
        % Low Pass
%          MedFiltR = 2.*(ceil(1/((f0Hz/2.5)*dtSec))-1)+1;
        % High Pass
%           MedFiltR = 2.*(ceil(1/((f0Hz*2.5)*dtSec))-1)+1;

        % Function
        Rad = medfilt1( Rad, MedFiltR, [], 2,'omitnan','truncate' );

%         Rad = Rad - MedRad;
        
        if isDisplay
        display( 'Median Subtraction Filter Done')
        toc
        display(' ')
        end
    end
    
    %----------------------------------------------------------------------
    % SVD First Principle Component Subtraction Filter
    if isSVDfilt
        if isDisplay
        display( 'Begin Singular Value Decomposition Filter')
        tic
        end
        
        Rad = SVDSfilter( Rad );
        
        if isDisplay
        display( 'Singular Value Decomposition Filter Done')
        toc
        display(' ')
        end
    end

    %----------------------------------------------------------------------    
    % Band-Pass Filter
    if isBandPass
        if isDisplay
        display( 'Begin Band-Pass Filter')
        tic
        end
        
        % Parameters
        % Build Two Octave Filter About Nominal Frequency
        fMin = f0Hz/2; % [Hz]   
        fMax = f0Hz*2; % [Hz]
%         fMin = f0Hz/2; % [Hz]   
%         fMax = f0Hz; % [Hz]
        
        % Function
        Rad = bpfilter( Rad, dtSec, fMin, fMax, 8 );
        
        if isDisplay
        display( 'Band-Pass Filter Done')
        toc
        display(' ')
        end
    end
            
    %----------------------------------------------------------------------
    % Time-Zero Correction
    if isTimeZero
        if isDisplay
        display( 'Begin Time Zero Correction')
        tic
        end
        
        % Parameters
        t = [0:dt:nsamp.*dt];
        merR = 1;   %Rank of MER window [ns]
        powMER = 3; % Power of MER operation
        % Function
        Rad = timeZero( Rad, t, dt, merR, powMER, offset);
        
        if isDisplay
        display( 'Time Zero Correction Done')
        toc
        display(' ')
        end
    end

    %----------------------------------------------------------------------
    % Exponential Time Dependant Gain
    if isExpGain
        if isDisplay
        display( 'Begin Power Gain')
        tic
        end
        
        % Parameters
        tpow = 2.25;%1.5; % Filter Order, 1 is Linear
        
        % Function
        Rad = gain(Rad, tpow );
        
        if isDisplay
        display( 'Power Gain Done')
        toc
        display(' ')
        end
    end
    
    % ---------------------------------------------------------------------
    % Automatic Gain Control
    if isACGGain
        if isDisplay
            tic
            display( 'Begin AGC')
        end
        
        % Parameters
%         R = f0/4.*dtSec; % Quarter Wavelength Convolution Filter
        R = ceil(nsamp/2);
        type = 2; % Trace Normalize: 0 = None, 1 = amplitude, 2 = RMSnorm

        Rad = AGCgain(Rad,R,type);
        
        if isDisplay
            display( 'AGC Done')
            toc
            display(' ')
        end
    end
    
    %----------------------------------------------------------------------
    % Clean Non-Existant Values
    if isCleanNaN
        if isDisplay
        tic
        display( 'Begin Clean NAN')
        end
        
        % Function
        [Rad,NaNno,datum] = cleanNaN(Rad);
        
        if isDisplay
        display( 'Clean NaN Done')
        toc
        cleanRecord = sprintf('Of %10.f Datum %10.f NaN Values Cleansed'...
            ,datum,NaNno);
        disp(cleanRecord)
        display(' ')
        end
    end
    
    %----------------------------------------------------------------------
    % Stacking Filter
    if isStak
        if isDisplay
        display( 'Begin Stack')
        tic
        end
        
        % Parameters
        StakFiltR = 10; % Filter Rank
        
        % Function
        Rad = StakR(Rad,StakFiltR);
        
        if isDisplay
        display( 'Stacking Done')
        toc
        display(' ')
        end
    end

end

