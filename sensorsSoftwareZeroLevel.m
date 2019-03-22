%% sensorsSoftwareZeroLevel.m
% This code unpacks the binary S&S data file and writes the raw traces and
% meta data to ncdf

% Tate Meehan - Oct. 2018
clear; close all; clc
addpath '/sonichome/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/TM'
isWrite = 1;

% Create Error Log
 f = fopen( 'ErrorLog.txt', 'w' );  
 fclose(f);
% Create Report Log
 f = fopen('ReportLog.txt','w');
%% From Data Directory Read and Sort SubDirectories
workingDirectory = pwd;
directory = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/PulseEKKO/500MHz';
folders = dir(directory);
folders(1:2) = []; % Remove Hidden Directories
datearray = cat(1,folders.date);
daymonthyear = (datearray(:,1:11));
num = datenum(daymonthyear,'dd-mmm-yyyy');
[~,sortIx] = sort(num);
folders = folders(sortIx);

% Loop Over each day of Acquisition
for ff = 16:length(folders)%1:length(folders)
    dataDir = [directory,'/',folders(ff).name];
    % This IS Robust for 2016 w/GPS files
    % Get .DT1 Files
    dt1filenames = dir(fullfile(dataDir,'*.DT1'));
    dt1files = struct2cell(dt1filenames);
    dt1files(2:end,:) = [];
    dt1files = cell2mat(dt1files');
    dt1bytes = [dt1filenames.bytes];
    % Remove Files Without Data from Query
    testbytes = dt1bytes == 0;
    if any(testbytes)
        rmDT1ix = find(testbytes);
        for ii = 1:length(any(testbytes))
            msg = [dataDir,'/',dt1files(rmDT1ix,:),' has 0 bytes.','\n'];
            warning(msg)
            WarnUser(msg,'ErrorLog.txt',workingDirectory)
        end
        dt1files(rmDT1ix,:) = [];
    end
    % Get .HD Files
    hdfilenames = dir(fullfile(dataDir,'*.HD'));
    hdfiles = struct2cell(hdfilenames);
    hdfiles(2:end,:) = [];
    hdfiles = cell2mat(hdfiles');
    
    % Construct Files
    files = dt1files(:,1:6);
    nFiles = size(files,1);
    
        if size(dt1files,1) ~= size(hdfiles,1)
            hdnum = str2num(hdfiles(:,5:6));
            nhd = hdnum+1;
            dt1num = str2num(dt1files(:,5:6));
            missingHD = dt1num(find(~ismember(dt1num,hdnum)));
            for ii = 1:length(missingHD)
                msg = ['Missing ', dataDir,'/',files(missingHD(ii),:),'.HD','\n'];
                warning(msg);
                WarnUser(msg,'ErrorLog.txt',workingDirectory)
            end
        else
        end
    
    for ii = 1 : nFiles
        tic
        %------------------------------------------------------------------
        % Multiplexed Channel Record
        filename = files(ii,:);
        
        % Read Data
        filepath = fullfile(dataDir,filename);
        [Rad,hdr1,trhd,dt,f0,~,dx,~] = readSensorsSoftwareData( filepath );
        % Append f0 to Trace Header
        trhd(9,:) = ones(1,size(trhd,2)).*f0;
        
        nChan = max(unique(trhd(23,:)));
        % Nominal Frequency GHz
        f0GHz = f0/1000;
        % No. Traces in Multiplexed Data
        [~, multiplexNtrcs] = size(Rad);
        
        % Install Transmitter and Receiver Sequencing and Geometry
        % 500 MHz Offsets
        if f0 == 500
            if nChan == 9
                % 9 Channel Offsets
                txGeo = [0, 0, 0, 1.33, 1.33, 1.33, 2.67, 2.67, 2.67]; % Tx Sequence & Absolute Position
                rxGeo = [4, 8, 12, 4, 8, 12, 4, 8, 12]; % Rx Sequence & Absolute Position
            elseif nChan == 6
                % 6 Channel Offsets
                txGeo = [1.33, 1.33, 1.33, 2.67, 2.67, 2.67]; % Tx Sequence & Absolute Position
                rxGeo = [4, 8, 12, 4, 8, 12]; % Rx Sequence & Absolute Position
            else
                display('Undefined Offset Array!')
            end
            % 1 GHz Offsets
        elseif f0 == 1000
            txGeo = [0, 0, 0, 0.67, 0.67, 0.67, 1.33, 1.33, 1.33]; % Tx Sequence & Absolute Position
            rxGeo = [2, 4, 6, 2, 4, 6, 2, 4, 6]; % Rx Sequence & Absolute Position
        end
        
        offsetArray = rxGeo - txGeo;            % Offset Array for CMP Gather
        % Append the Offset Array for NetCDF Export
        channelArray = trhd(23,:);
        offsetAppend = offsetArray(channelArray);
        trhd(22,:) = offsetAppend;
        
        % Attach trace header to data for NetCDF Export
        DATA = [trhd;Rad];
        
        %         GPR = struct('metaData',{hdr1},'traceHeader',trhd,'traces',Rad,...
        %             'offsets',offsetArray,'frequency',f0,'timeSample',dt,'spaceSample',dx);
        
        %% Write to NetCDF
        if isWrite
            % Create the output Directory if it does not exist
            if exist([workingDirectory,'/',folders(ff).name],'dir') ~= 7
                mkdir(workingDirectory,folders(ff).name)
            end
            % Change to outdir for File Write
            cd([workingDirectory,'/',folders(ff).name])
            
            % Name the NetCDF.nc file
            ncdfName = [folders(ff).name,'-',filename,'.nc'];
            
            fprintf(['Writing ',ncdfName, '\n'])
            
            numrow = size(DATA,1);
            
            numcol = size(DATA,2);
            
            % Create the NetCDF-4 file with Dimensions of DATA
            netcdf.setDefaultFormat('FORMAT_NETCDF4') ; 
            
            ncid = netcdf.create(ncdfName,'CLOBBER');
            
            dimidrow = netcdf.defDim(ncid,'rows',numrow);
            
            dimidcol = netcdf.defDim(ncid,'length',numcol);
            
            % Define and put Variable DATA
            varid = netcdf.defVar(ncid,'DATA','NC_DOUBLE',[dimidrow dimidcol]);
            
            % Increase NetCDF Compression
            netcdf.defVarDeflate(ncid,varid,true,true,5);
            
            netcdf.endDef(ncid);
            
            netcdf.putVar(ncid,varid,DATA);
                        
            % Close the .nc file after writing
            netcdf.close(ncid);
            % Test for Proper Writing
            %     ncid2 = netcdf.open(ncdfName,'NC_NOWRITE');
            %
            %     data_copy = netcdf.getVar(ncid2,0);
            %
            %     if isequal(DATA,data_copy)
            %
            %         disp('Data match');
            %
            %     else
            %
            %         disp('Data mis-match');
            %     end
            
            % Return to the working directory for next iteration
            cd(workingDirectory)
            fprintf(['Done writing ',folders(ff).name,'-',files(ii,:),'\n'])
            reportMessage = ['Wrote ',folders(ff).name,'-',filename,'.nc','\n'];
            Report2User(reportMessage,'ReportLog.txt',workingDirectory)
        end        
        toc
        display(' ')
    end
end