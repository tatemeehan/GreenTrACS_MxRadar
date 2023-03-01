function [GPR] = readMXnc(GPR)
%     ii = GPR.MD.ii;
for ii = 1 : GPR.MD.nFiles
    %------------------------------------------------------------------
    % Multiplexed Channel Record
%     filename = GPR.MD.fileNames(ii).name;
%     filename = GPR.MD.Dir.name;
    filename = GPR.MD.fileNames;
    filepath = fullfile(GPR.MD.dataDir,filename);
    % Read netCDF data file
    disp(' ')
    fprintf('Reading .nc File \n')
    tic
    ncRad = ncread(filepath,'DATA');
    GPR.D.trhd{ii} = ncRad(1:27,:);
    GPR.D.MxRadar{ii} = ncRad(28:end,:);
    clear('ncRad');
    % Load Geometry
    GPR.Geometry.Chan{ii} = unique(GPR.D.trhd{ii}(3,:));
    GPR.Geometry.nChan{ii} = numel(unique(GPR.D.trhd{ii}(3,:)));
    GPR.Geometry.offset{ii} = GPR.D.trhd{ii}(4,1:GPR.Geometry.nChan{ii});
    GPR.Geometry.midpointx{ii} = zeros(1,GPR.Geometry.nChan{ii});
    GPR.Geometry.midpointy{ii} = GPR.D.trhd{ii}(5,1:GPR.Geometry.nChan{ii});
%     GPR.Geometry.tx{ii} = GPR.D.trhd{ii}(7,1:GPR.Geometry.nChan{ii});
%     GPR.Geometry.ty{ii} = GPR.D.trhd{ii}(8,1:GPR.Geometry.nChan{ii});
%     GPR.Geometry.rx{ii} = GPR.D.trhd{ii}(9,1:GPR.Geometry.nChan{ii});
%     GPR.Geometry.ry{ii} = GPR.D.trhd{ii}(10,1:GPR.Geometry.nChan{ii});
%     GPR.Geometry.GPSxyz{ii} = GPR.D.trhd{ii}(35:37,1);
%     GPR.Geometry.GPSdelta{ii} = GPR.D.trhd{ii}(38:40,1:GPR.Geometry.nChan{ii});
    
    % Load Geolocation
    GPR.Geolocation.X{ii} = GPR.D.trhd{ii}(10,1:GPR.Geometry.nChan{ii}:end);
    GPR.Geolocation.Y{ii} = GPR.D.trhd{ii}(11,1:GPR.Geometry.nChan{ii}:end);
    GPR.Geolocation.Z{ii} = GPR.D.trhd{ii}(12,1:GPR.Geometry.nChan{ii}:end);
    GPR.Geolocation.Longitude{ii} = GPR.D.trhd{ii}(22,1:GPR.Geometry.nChan{ii}:end);
    GPR.Geolocation.Latitude{ii} = GPR.D.trhd{ii}(23,1:GPR.Geometry.nChan{ii}:end);
    GPR.Geolocation.Elevation{ii} = GPR.D.trhd{ii}(24,1:GPR.Geometry.nChan{ii}:end);
    GPR.Geolocation.Distance{ii} = GPR.D.trhd{ii}(13,1:GPR.Geometry.nChan{ii}:end);
    GPR.Geolocation.Speed{ii} = GPR.D.trhd{ii}(14,1:GPR.Geometry.nChan{ii}:end);
    GPR.Geolocation.Heading{ii} = GPR.D.trhd{ii}(15,1:GPR.Geometry.nChan{ii}:end);
    GPR.Geolocation.Slope{ii} = GPR.D.trhd{ii}(16,1:GPR.Geometry.nChan{ii}:end);
    GPR.Geolocation.Geoid{ii} = GPR.D.trhd{ii}(24,1:GPR.Geometry.nChan{ii}:end)-GPR.D.trhd{ii}(12,1:GPR.Geometry.nChan{ii}:end);
%     GPR.Geolocation.UTMzone{ii} = GPR.D.trhd{ii}(34,1:GPR.Geometry.nChan{ii}:end);

    % Load GPR Parameters
    GPR.D.f0{ii} = (GPR.D.trhd{ii}(6,1)); % [MHz]
    GPR.D.f0GHz{ii} = GPR.D.f0{ii}/1000; %[GHz]
    GPR.D.dt{ii} = GPR.D.trhd{ii}(7,1); % [ns]
    GPR.D.dx{ii} = mean(diff(GPR.Geolocation.Distance{ii})); % [m]
    GPR.D.TimeAxis{ii} = [0:GPR.D.dt{ii}:(GPR.D.trhd{ii}(8,1)-1).*GPR.D.dt{ii}]';
    
    % Load Meta Data
    d = datetime(GPR.D.trhd{ii}(9,:),'ConvertFrom','posixtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
    GPR.MD.DateTime{ii} = d;
    disp(' ')
    fprintf('Read .nc File \n')
    tic
end
end

