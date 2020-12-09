close all
% clear
clc

verbose=0;

%% Get gage location (Latitude 35°46'25", Longitude 113°21'46")
% Source: https://waterdata.usgs.gov/nwis/inventory?agency_code=USGS&site_no=09404200
gage=[dms2degrees([35,46,25]), -dms2degrees([113,21,46])];

% Source: https://waterdata.usgs.gov/monitoring-location/09404200/#parameterCode=00060
gage_precise=[35.7735994, -113.363544];

% Convert (lat,lon) to (Easting,Northing)
dczone = utmzone(gage(1),gage(2));
utmstruct = defaultm('utm');
utmstruct.zone = dczone;
utmstruct.geoid = referenceEllipsoid('Geodetic Reference System 1980');
utmstruct = defaultm(utmstruct);
[E_gage,N_gage] = mfwdtran(utmstruct,gage(1),gage(2));

if verbose
    figure
    hold on
    plot(E_gage,N_gage,'*')

    figure
    geoshow(gage(1),gage(2),'DisplayType','point')
end

%% Get location of the vessel from data survey files (XTF)
% Version 1: considering the location provided by POSRAW
if verbose
    t=0; 
    name4legend={}; % just name of which xtf for displaying in the legends
    figure(1)
    hold on
    grid on

    figure(2)
    hold on
    grid on
end

trajectory_POSRAW={};
k=1;

for iRun=['A','B','C','D']
    cd(strcat('/Users/thanhhuynguyen/Storage/DATA/02_DiamondCreek/mat_truepix/Run_',iRun,'/'))
    list=dir('*.mat');

    for iFile=1:length(list)
        load(list(iFile).name)
        
        % Get lat and lon from posraw packet
        lon=[];
        lat=[];
        for i=1:length(Ping)
            if Ping(i).HeaderType==107
                lon=[lon Ping(i).xtf_posraw_navigation.RawXcoordinate];
                lat=[lat Ping(i).xtf_posraw_navigation.RawYcoordinate];
            end
        end

        % Convert (lat,lon) to (Easting,Northing)
        dczone = utmzone(mean(lat,'omitnan'),mean(lon,'omitnan'));
        utmstruct = defaultm('utm');
        utmstruct.zone = dczone;
        utmstruct.geoid = referenceEllipsoid('Geodetic Reference System 1980');
        utmstruct = defaultm(utmstruct);
        [E,N] = mfwdtran(utmstruct,lat,lon);
        
        if verbose
            figure(1)
            plot(E,N,'-','color',rand(1,3))

            figure(2)
            geoshow(lat,lon,'DisplayType','multipoint')

            t=t+1;
            name4legend(t,:)={strcat(iRun,' - ',list(iFile).name)};
        end
        
        trajectory_POSRAW(k,1)={strcat(iRun,num2str(iFile))};
        trajectory_POSRAW(k,2)={lat};
        trajectory_POSRAW(k,3)={lon};
        k=k+1;
    end
end
if verbose
    figure(1)
    hold on
    plot(E_gage,N_gage,'*')
    legend(name4legend,'Gage')

    figure(2)
    hold on
    geoshow(gage(1),gage(2),'DisplayType','point')
    legend(name4legend,'Gage')
end

%% Version 2: considering the location provided by XTFBathyHeader
if verbose
    t=0; 
    name4legend={}; % just name of which xtf for displaying in the legends
    figure(1)
    hold on
    grid on

    figure(2)
    hold on
    grid on
end

trajectory_bathy_header={};
k=1;
for iRun=['A','B','C','D']
    cd(strcat('/Users/thanhhuynguyen/Storage/DATA/02_DiamondCreek/mat_truepix/Run_',iRun,'/'))
    list=dir('*.mat');

    for iFile=1:length(list)
        load(list(iFile).name)

        % Get lat and lon from xtfbathy_header packet
        lon=[];
        lat=[];
        for i=1:length(Ping)
            if Ping(i).HeaderType==67
                lon=[lon Ping(i).xtfbathy_header.ShipXcoordinate];
                lat=[lat Ping(i).xtfbathy_header.ShipYcoordinate];
            end
        end

        % Convert (lat,lon) to (Easting,Northing)
        dczone = utmzone(mean(lat,'omitnan'),mean(lon,'omitnan'));
        utmstruct = defaultm('utm');
        utmstruct.zone = dczone;
        utmstruct.geoid = referenceEllipsoid('Geodetic Reference System 1980');
        utmstruct = defaultm(utmstruct);
        [E,N] = mfwdtran(utmstruct,lat,lon);

        if verbose
            figure(1)
            plot(E,N,'-','color',rand(1,3))

            figure(2)
            geoshow(lat,lon,'DisplayType','multipoint')

            t=t+1;
            name4legend(t,:)={strcat(iRun,' - ',list(iFile).name)};
        end
        
        trajectory_bathy_header(k,1)={strcat(iRun,num2str(iFile))};
        trajectory_bathy_header(k,2)={lat};
        trajectory_bathy_header(k,3)={lon};
        k=k+1;
    end
end
if verbose
    figure(1)
    hold on
    plot(E_gage,N_gage,'*')
    legend(name4legend,'Gage')

    figure(2)
    hold on
    geoshow(gage(1),gage(2),'DisplayType','point')
    legend(name4legend,'Gage')
end

%% Set a look-up distance from gage location 
% + Extract the packets in the region-of-interest
% + Save the extracted results

tol=100; % allowed tolerance 
E_gage_tol=E_gage+[-tol tol]; % look-up distance from gage
N_gage_tol=N_gage+[-tol tol];

for iRun=['A','B','C','D']
    cd(strcat('/Users/thanhhuynguyen/Storage/DATA/02_DiamondCreek/mat_truepix/Run_',iRun,'/'))
    list=dir('*.mat');

    for iFile=1:length(list)
        load(list(iFile).name)
    
        tbPing107=[]; % table of Pings 107
        tbPingNearGage=[]; % table of Pings near gage

        lon=[];
        lat=[];
        for i=1:length(Ping)
            if Ping(i).HeaderType==107
                lon=[lon Ping(i).xtf_posraw_navigation.RawXcoordinate];
                lat=[lat Ping(i).xtf_posraw_navigation.RawYcoordinate];
                
                tbPing107=[tbPing107 i];
            end
        end

        % Convert (lat,lon) to (Easting,Northing)
        dczone = utmzone(mean(lat,'omitnan'),mean(lon,'omitnan'));
        utmstruct = defaultm('utm');
        utmstruct.zone = dczone;
        utmstruct.geoid = referenceEllipsoid('Geodetic Reference System 1980');
        utmstruct = defaultm(utmstruct);
        [E,N] = mfwdtran(utmstruct,lat,lon);

        
        for i=1:length(tbPing107)
            if (E(i)>E_gage_tol(1) && E(i)<E_gage_tol(2)) &&...
                    (N(i)>N_gage_tol(1) && N(i)<N_gage_tol(2))
                disp(['Run ', iRun, ' ping=',num2str(tbPing107(i))])
                
                % Save the Ping 107 in the ROI and 3 subsequent Pings after it
                tbPingNearGage=[tbPingNearGage tbPing107(i):(tbPing107(i)+3)];   
            end
        end
        
        % Removing other pings and Saving the extraction result
        if ~isempty(tbPingNearGage)
            Ping_new=Ping;
            tbPingsToRemove=setdiff(1:length(Ping), tbPingNearGage);
            Ping_new(tbPingsToRemove)=[];

            survey_name=strsplit(list(iFile).name,'.xtf.mat');
            save(strcat(survey_name{1},'_nearGage.mat'),'Ping_new','xtffile_header')
        end
    end
end

%% Write shapefiles (added on 21/10/2020)
dir_data='/Users/thanhhuynguyen/Storage/DATA/02_DiamondCreek/Trajectory_SHP/';
for i=1:length(trajectory_POSRAW)
    lat=trajectory_POSRAW{i,2};
    lon=trajectory_POSRAW{i,3};
    s = geoshape();
    s.Latitude = lat;
    s.Longitude = lon;
    
    shapewrite(s,strcat(dir_data,'POSRAW_',trajectory_POSRAW{i,1},'.shp'))
end

for i=1:length(trajectory_bathy_header)
    lat=trajectory_bathy_header{i,2};
    lon=trajectory_bathy_header{i,3};
    s = geoshape();
    s.Latitude = lat;
    s.Longitude = lon;
    
    shapewrite(s,strcat(dir_data,'bathy_header_',trajectory_bathy_header{i,1},'.shp'))
end

s = geopoint();
s.Latitude = gage(1);
s.Longitude = gage(2);

shapewrite(s,strcat(dir_data,'GAGE_location','.shp'))

s = geopoint();
s.Latitude = gage_precise(1);
s.Longitude = gage_precise(2);

shapewrite(s,strcat(dir_data,'GAGE_location_precise_','.shp'))


