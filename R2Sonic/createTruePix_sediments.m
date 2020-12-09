% Created by TH Nguyen, ULaval
% Last modified on: 2019/11/07
%
%% Notes
% Each .mat file corresponds to one XTF file.
% A .mat file contains an "xtf_fileheader" and a "Ping". Both are struct
% variables containing multiple fields.
% The structure of both of them are exactly in the order of decoding, and
% the name of fields are set exactly like in the XTF documentation. A field
% can also be another struct, which contains other children fields.
%
% For example:
% To get the value of PacketName in the 3rd ping (a Truepix ping) use:
% >> Ping(3).truepix.PacketName
%
% Output: Truepix data converted into a 3D point cloud (with intensity)
%         Outputs are stocked in <dirResult>/mat/<filename>_pts.mat
% 
%% Load .mat file
close all
clc

verbose=0;
currentFolder = pwd;
% data_dir = '/Users/thanhhuynguyen/Storage/r2SonicChallenge19/XTF_truepix/mat/'
% dirResult = '/Users/thanhhuynguyen/Storage/r2SonicChallenge19/XTF_truepix2/'
% data_dir = '/Users/thanhhuynguyen/Storage/r2SonicChallenge19/Eastcom_20190529 100-200-400kHz_stationary/'
% dirResult = '/Users/thanhhuynguyen/Storage/r2SonicChallenge19/Eastcom_20190529 100-200-400kHz_stationary/'
% data_dir = '/Users/thanhhuynguyen/Storage/r2SonicChallenge19/DC_20190226/mat/'

data_dir = '/Users/thanhhuynguyen/Storage/DATA/02_DiamondCreek/mat_truepix/Run A/'
dirResult = data_dir

cd(data_dir)
list = dir('*.xtf.mat');

cd(currentFolder)
for idxFile=1:size(list,1)
    
    fileName = list(idxFile).name
    load(strcat(data_dir,fileName))
    
    %% Read POSRAW_NAVIGATION packets (HeaderType=107)
    time_nav=[];
    posraw_nav=[];
    for iPing=1:length(Ping)
        if Ping(iPing).HeaderType==107
            time_nav=[time_nav;
                (Ping(iPing).xtf_posraw_navigation.Hour*3600+...
                Ping(iPing).xtf_posraw_navigation.Minute*60+...
                Ping(iPing).xtf_posraw_navigation.Second)*1e3+...
                Ping(iPing).xtf_posraw_navigation.MicroSeconds/10];
            
            posraw_nav=[posraw_nav;
                Ping(iPing).xtf_posraw_navigation.RawXcoordinate...
                Ping(iPing).xtf_posraw_navigation.RawYcoordinate...
                Ping(iPing).xtf_posraw_navigation.RawAltitude];
        end
    end
    
    if (verbose)
        figure
        plot3(posraw_nav(:,1), posraw_nav(:,2), posraw_nav(:,3))
        grid on
        xlabel('Longitude [°]')
        ylabel('Lattitude [°]')
        zlabel('Altitude [m]')
    end
    
    
    %% Read ATTITUDE_DATA packets (HeaderType=3)
    % Note: the data provided by ATTITUDE_DATA packets are to be used as
    % the attitude of the sonar at the reception of a ping, computed
    % according to the TimeTag.
    
    timetagRx=[];
    rollRx=[];
    pitchRx=[];
    heaveRx=[];
    yawRx=[];
    headingRx=[];
    
    for iPing=1:length(Ping)
        if Ping(iPing).HeaderType==3
            timetagRx=[timetagRx; Ping(iPing).attitudedata.TimeTag];
            
            rollRx=[rollRx; Ping(iPing).attitudedata.Roll];
            pitchRx=[pitchRx; Ping(iPing).attitudedata.Pitch];
            yawRx=[yawRx; Ping(iPing).attitudedata.Yaw];
            heaveRx=[heaveRx; Ping(iPing).attitudedata.Heave];
            headingRx=[headingRx; Ping(iPing).attitudedata.Heading];
        end
    end
    
    if (verbose)
        figure
        subplot(221), plot(timetagRx, rollRx), xlabel('TimeTag'), ylabel('Roll [°]'), grid on
        subplot(222), plot(timetagRx, pitchRx), xlabel('TimeTag'), ylabel('Pitch [°]'), grid on
        subplot(223), plot(timetagRx, heaveRx), xlabel('TimeTag'), ylabel('Heave [m]'), grid on
        subplot(224), plot(timetagRx, headingRx+yawRx), xlabel('TimeTag'), ylabel('Heading+Yaw [°]'), grid on
    end
    
    %% Read attitude data from the truepix data packet
    % Note: these data are attitude of the sonar at the ping emission
    timetagTx=[];
    rollTx=[];
    pitchTx=[];
    yawTx=[];
    heaveTx=[];
    headingTx=[];
    
    nbPing67=0;  % nb of ping 67 (i.e. truepix packets)
    idxPing67=[]; % array of indices of ping 67 (i.e. truepix packets)
    
    for iPing=1:length(Ping)
        if Ping(iPing).HeaderType==67
            
            nbPing67  = nbPing67+1;
            idxPing67 = [idxPing67, iPing];
            
            timetagTx=[timetagTx; Ping(iPing).xtfbathy_header.AttitudeTimeTag];
            
            rollTx=[rollTx; Ping(iPing).xtfbathy_header.SensorRoll];
            pitchTx=[pitchTx; Ping(iPing).xtfbathy_header.SensorPitch];
            yawTx=[yawTx; Ping(iPing).xtfbathy_header.Yaw];
            heaveTx=[heaveTx; Ping(iPing).xtfbathy_header.Heave];
            headingTx=[headingTx; Ping(iPing).xtfbathy_header.SensorHeading];
        end
    end
    
    if (verbose)
        figure
        subplot(221), plot(timetagTx, rollTx), xlabel('TimeTag'), ylabel('Roll [°]'), grid on
        subplot(222), plot(timetagTx, pitchTx), xlabel('TimeTag'), ylabel('Pitch [°]'), grid on
        subplot(223), plot(timetagTx, heaveTx), xlabel('TimeTag'), ylabel('Heave [m]'), grid on
        subplot(224), plot(timetagTx, headingTx+yawTx), xlabel('TimeTag'), ylabel('Heading [°]'), grid on
    end
    
    %% Interpolate/Calculate actual location and attitude for each ping
    % Note:
    %    For location, use the POSRAW_NAVIGATION packets as a function of
    %    time to interpolate the location of the sonar at emission moment
    
    x=interp1(time_nav,posraw_nav(:,1),timetagTx);
    y=interp1(time_nav,posraw_nav(:,2),timetagTx);
    z=interp1(time_nav,posraw_nav(:,3),timetagTx);
    posraw_nav67=[x y z];
    
    % Note:
    %    For roll, we use the value at the reception of a ping. For pitch,
    %    we use the value at ping emission. For yaw heave and heading, we
    %    use the average value of the sonar's yaw at emission and reception
    
    pitch67=pitchTx;
    
    roll67=zeros(1,nbPing67);
    yaw67=zeros(1,nbPing67);
    heave67=zeros(1,nbPing67);
    heading67=zeros(1,nbPing67);
    
    for i=1:nbPing67
        iPing=idxPing67(i);
        
        c=double(Ping(iPing).truepix.H0_SoundSpeed);
        r=double(Ping(iPing).truepix.H0_RxRange);
        delta_t=2*r/c;
        
        roll67(i)   =interp1(timetagRx,rollRx,timetagTx(i)+delta_t*1000);
        yaw67(i)    =0.5*(interp1(timetagRx,yawRx,timetagTx(i)+delta_t*1000) + yawTx(i));
        heave67(i)  =0.5*(interp1(timetagRx,heaveRx,timetagTx(i)+delta_t*1000) + heaveTx(i));
        heading67(i)=0.5*(interp1(timetagRx,headingRx,timetagTx(i)+delta_t*1000) + headingTx(i));
    end
    
    
    %% Store magnitude and angle data (Truepix data, type D1)
    
    port_mag=[]; % Magnitude data from port side
    stbd_mag=[]; % Magnitude data from starboard side
    port_ang=[]; % Angle data from port side (multiplied with scale factor)
    stbd_ang=[]; % Angle data from starboard side (multiplied with scale factor)
    for iPing=1:length(Ping)
        if Ping(iPing).HeaderType==67
            port_mag=[port_mag; Ping(iPing).truepix.D_all_packets(:,1)'];
            stbd_mag=[stbd_mag; Ping(iPing).truepix.D_all_packets(:,3)'];
            
            angle_scaling_factor=Ping(iPing).truepix.D(1).AngleScalingFactor;
            port_ang=[port_ang; Ping(iPing).truepix.D_all_packets(:,2)'.*angle_scaling_factor];
            stbd_ang=[stbd_ang; Ping(iPing).truepix.D_all_packets(:,4)'.*angle_scaling_factor];
        end
    end
    
    if (verbose)
        figure
        subplot(221), plot(Ping(3).truepix.D_all_packets(:,1)), grid on, xlabel('Samples'), ylabel('Port Mag. [\muPa]')
        subplot(223), plot(Ping(3).truepix.D_all_packets(:,2).*angle_scaling_factor), grid on, xlabel('Samples'), ylabel('Port Angle [rad]')
        subplot(222), plot(Ping(3).truepix.D_all_packets(:,3)), grid on, xlabel('Samples'), ylabel('Stbd Mag. [\muPa]')
        subplot(224), plot(Ping(3).truepix.D_all_packets(:,4).*angle_scaling_factor), grid on, xlabel('Samples'), ylabel('Stbd Angle [rad]')
    end
    
    %% Requantization (16-bit data to 8-bit data) -> For display purpose
    C=20;
    port_mag_req=C*log10(1+port_mag(:,:)*2^8/2^16);
    stbd_mag_req=C*log10(1+stbd_mag(:,:)*2^8/2^16);
    
    if (verbose)
        figure
        subplot(121), imagesc(fliplr(port_mag_req)), xlabel('Port side Magnitude')
        subplot(122), imagesc(stbd_mag_req), xlabel('Starboard side Magnitude')
        
        figure
        subplot(121), imagesc(fliplr(port_ang)), xlabel('Port side Angle')
        subplot(122), imagesc(stbd_ang), xlabel('Starboard side Angle')
    end
    
    
    %% Applying TVG
    tvg=zeros(nbPing67,double(Ping(3).truepix.D(1).TotalSamples));
    for i = 1:nbPing67
        iPing=idxPing67(i)
        
        fs=Ping(iPing).truepix.H0_RxSampleRate;
        c=Ping(iPing).truepix.H0_SoundSpeed;
        range=[1:double(Ping(iPing).truepix.D(1).TotalSamples)]*double(c/2/fs);
        
        alpha=Ping(iPing).truepix.H0_RxAbsorption; % Absorption
        Sp=Ping(iPing).truepix.H0_RxSpreading; % Spreading rate
        G=2*Ping(iPing).truepix.H0_RxGain; % from data format doc -> multiply by 2 for dB
        tvg(i,:)=2*range*alpha/1000 + Sp*log(range) + G;
    end
    
    port_mag_tvg=zeros(size(port_mag));
    for i = 1:nbPing67
        port_mag_tvg(i,:)=port_mag(i,:).*10.^(tvg(i,:)/10);
    end
    port_mag_tvg_req=C*log(1+port_mag_tvg(:,:)*2^8/2^16);
    
    %     figure
    %     subplot(121), imagesc(port_mag_req(:,475:6290))
    %     colormap(gray)
    %     subplot(122), imagesc(port_mag_tvg_req(:,475:6290),[30 160])
    %     colormap(gray)
    %
    %% Convert ship trajectory from (long, lat) into UTM coordinate system
    lon=posraw_nav67(:,1);
    lat=posraw_nav67(:,2);
%     Z=posraw_nav67(:,3);
    
    dczone = utmzone(mean(lat,'omitnan'),mean(lon,'omitnan'));
    
    utmstruct = defaultm('utm');
    utmstruct.zone = dczone;
    
    if strcmp(xtffile_header.SpheriodType(1:8),'GRS 1980')
        utmstruct.geoid = referenceEllipsoid('Geodetic Reference System 1980');
    elseif strcmp(xtffile_header.SpheriodType(1:8),'WGS 1984')
        utmstruct.geoid = referenceEllipsoid('wgs84');
    end
    
    utmstruct = defaultm(utmstruct);
    
%     [E,N] = mfwdtran(utmstruct,lat,lon);
    [E,N,Z] = mfwdtran(utmstruct,lat,lon,posraw_nav67(:,3));
    
    % figure
    % plot3(E,N,Z,'*'), grid on
    
    x0 = E; % Easting (in meters)
    y0 = N; % Northing (in meters)
    z0 = Z; % Ship altitude
    
    ship_trajectory = [x0 y0 z0];
    
    if (verbose)
        figure
        plot3(x0, y0, z0, 'rv')
        hold on, grid on
        for i=1:5:nbPing67
            text(x0(i),y0(i), z0(i), num2str(i))
        end
        view(2)
    end
    
    %% Convert the truepix data (magnitude, angle, range) into 3D point cloud
    
    kappa = heading67+yaw67; % Heading angle + Yaw
    phi = pitch67;
    omega = roll67;
    
    % lever arm translation
    a_bI=[xtffile_header.ChanInfo(1).OffsetX,...
        xtffile_header.ChanInfo(1).OffsetY,...
        xtffile_header.ChanInfo(1).OffsetZ];
    
    % boresight angles
%     ang_boresight=[xtffile_header.ChanInfo(1).OffsetRoll,...
%         xtffile_header.ChanInfo(1).OffsetPitch,...
%         xtffile_header.ChanInfo(1).OffsetYaw];
    ang_boresight=[0.910, 5.810, 2.980];
    r_bS=rotz(ang_boresight(3))*roty(ang_boresight(2))*rotx(ang_boresight(1));

    points_all=[];
    for i = 1:nbPing67-1
        iPing=idxPing67(i)
        
        fs=Ping(iPing).truepix.H0_RxSampleRate;
        c=Ping(iPing).truepix.H0_SoundSpeed;
        
        iSample=1:double(Ping(iPing).truepix.D(1).TotalSamples);
        range=iSample*double(c/2/fs);
        
        % samples from port side (in the sensor's frame of reference)
        p1=[zeros(length(iSample),1),...
            range(iSample)'.*sin(port_ang(i,iSample)'),...
            -range(iSample)'.*cos(port_ang(i,iSample)')];
        
        % samples from stbd side (in the sensor's frame of reference)
        p2=[zeros(length(iSample),1),...
            range(iSample)'.*sin(stbd_ang(i,iSample)'),...
            -range(iSample)'.*cos(stbd_ang(i,iSample)')];
        
        p_bS=[p1; p2]; % Each row is [0, Rcos(theta), -Rsin(theta)]
        
        P_prime = rotz(-kappa(i)+90)*roty(-phi(i))*rotx(-omega(i))*(r_bS*p_bS'-a_bI');
        
        points_all=[points_all; ship_trajectory(i,:) + P_prime'];
    end
    
    if (verbose)
        figure
        pcshow(points_all)
    end
    
    %%
    mag_thres=200;
    points_thres=[];
    points_thres_with_int=[];
    for i = 1:nbPing67-1
        iPing=idxPing67(i)
        
        fs=Ping(iPing).truepix.H0_RxSampleRate;
        c=Ping(iPing).truepix.H0_SoundSpeed;
        
        iSample=1:double(Ping(iPing).truepix.D(1).TotalSamples);
        range=iSample*double(c/2/fs);
        
        % select only the samples with significant magnitude (i.e. > mag_thres)
        iSample_port=find(port_mag(i,:)>mag_thres);
        iSample_stbd=find(stbd_mag(i,:)>mag_thres);
        
        p1=[zeros(length(iSample_port),1),...
            range(iSample_port)'.*sin(port_ang(i,iSample_port)'),...
            -range(iSample_port)'.*cos(port_ang(i,iSample_port)')];
        
        p2=[zeros(length(iSample_stbd),1),...
            range(iSample_stbd)'.*sin(stbd_ang(i,iSample_stbd)'),...
            -range(iSample_stbd)'.*cos(stbd_ang(i,iSample_stbd)')];
        
        p_bS=[p1; p2];
        
        P_prime = rotz(-kappa(i)+90)*roty(-phi(i))*rotx(-omega(i))*(r_bS*p_bS'-a_bI');
        
        points_thres=[points_thres; ship_trajectory(i,:) + P_prime'];
        
        points_thres_with_int=[points_thres_with_int;
            [ship_trajectory(i,:)+P_prime'...
            [port_mag_req(i,iSample_port) stbd_mag_req(i,iSample_stbd)]']];
        
    end
    
    if (verbose)
        figure
        pcshow(points_thres)
        
        figure
        pcshow(points_thres_with_int(:,1:3),points_thres_with_int(:,4))
    end
    
    %% Result: point cloud with only significant magnitude points
    if (verbose)
        figure(idxFile)
        pcshow(points_thres_with_int(:,1:3),points_thres_with_int(:,4))
        hold on
        plot3(E,N,Z,'r>'), grid on
        legend('3D points','Sonar position')
        xlabel('Easting (UTM)')
        ylabel('Northing (UTM)')
        zlabel('Altitude')
    end
    
    s=strsplit(fileName,'.');
%     save(strcat(dirResult,s{1},'_pts.mat'),'points_all','points_thres_with_int','ship_trajectory')
    
    
end


%% Display results (all truepix point clouds together)
% if (verbose)
    for idxFile=1:size(list,1)
        
        fileName = list(idxFile).name
        s=strsplit(fileName,'.');
        load(strcat(dirResult,s{1},'_pts.mat'),'points_thres_with_int','ship_trajectory')
        
        figure(6)
        color=rand(3,1);
%         pcshow(points_thres_with_int(:,1:3))
        pcshow(points_thres_with_int(:,1:3),points_thres_with_int(:,4))
%         plot3(points_thres_with_int(:,1),points_thres_with_int(:,2),points_thres_with_int(:,3),'.','Color',color)
        hold on
        plot3(ship_trajectory(:,1),ship_trajectory(:,2),ship_trajectory(:,3),'>','Color',color) 
        grid on
        legend('3D points','Sonar position')
        xlabel('Easting (UTM)')
        ylabel('Northing (UTM)')
        zlabel('Altitude')
    end
% end
%%

nbPing67=0;  % nb of ping 67 (i.e. truepix packets)
idxPing67=[]; % array of indices of ping 67 (i.e. truepix packets)

for iPing=1:length(Ping)
    if Ping(iPing).HeaderType==67
        
        nbPing67  = nbPing67+1;
        idxPing67 = [idxPing67, iPing];
    end
end

%%
idx170=[];
idx330=[];
idx450=[];
idx700=[];

for i=1:nbPing67
    iPing=idxPing67(i);

    if(Ping(iPing).truepix.H0_Frequency==170000)
        idx170=[idx170 i];
    elseif(Ping(iPing).truepix.H0_Frequency==330000)
        idx330=[idx330 i];
    elseif(Ping(iPing).truepix.H0_Frequency==450000)
        idx450=[idx450 i];
    else
        idx700=[idx700 i];
    end
end


for i = 1:10%nbPing67
    iPing=idxPing67(i);
    disp(Ping(iPing).truepix.H0_Frequency)
    fs=Ping(iPing).truepix.H0_RxSampleRate;
    c=Ping(iPing).truepix.H0_SoundSpeed;
    
%     disp(Ping(iPing).truepix.H0_TxBeamwidthHoriz)
%     disp(Ping(iPing).truepix.H0_TxBeamwidthVert)
    double(Ping(iPing).truepix.D(1).TotalSamples)
%     range=[1:double(Ping(iPing).truepix.D(1).TotalSamples)]*double(c/2/fs);
    
end
%%

figure
subplot(121), imagesc(fliplr(port_mag_req(idx700,:))), xlabel('Port side Magnitude')
subplot(122), imagesc(stbd_mag_req(idx700,:)), xlabel('Starboard side Magnitude')

figure
subplot(121), imagesc(fliplr(port_mag(idx700,:))), xlabel('Port side Magnitude')
subplot(122), imagesc(stbd_mag(idx700,:)), xlabel('Starboard side Magnitude')

figure
subplot(221), imagesc(fliplr(port_mag_req(idx170,:))), xlabel('Port side Magnitude (170 kHz)')
subplot(222), imagesc(fliplr(port_mag_req(idx330,:))), xlabel('Port side Magnitude (330 kHz)')
subplot(223), imagesc(fliplr(port_mag_req(idx450,:))), xlabel('Port side Magnitude (450 kHz)')
subplot(224), imagesc(fliplr(port_mag_req(idx700,:))), xlabel('Port side Magnitude (700 kHz)')

figure
subplot(221), imagesc(fliplr(port_ang(idx170,:))), xlabel('Port side Magnitude (170 kHz)')
subplot(222), imagesc(fliplr(port_ang(idx330,:))), xlabel('Port side Magnitude (330 kHz)')
subplot(223), imagesc(fliplr(port_ang(idx450,:))), xlabel('Port side Magnitude (450 kHz)')
subplot(224), imagesc(fliplr(port_ang(idx700,:))), xlabel('Port side Magnitude (700 kHz)')


%% Adding ship trajectory, using UTM coordinate system

lat=posraw_nav67(:,2);
lon=posraw_nav67(:,1);
Z=posraw_nav67(:,3);

dczone = utmzone(mean(lat,'omitnan'),mean(lon,'omitnan'));

utmstruct = defaultm('utm');
utmstruct.zone = dczone;
utmstruct.geoid = wgs84Ellipsoid;
utmstruct = defaultm(utmstruct);

[E,N] = mfwdtran(utmstruct,lat,lon);

% figure
% plot3(E,N,Z,'*'), grid on

x0 = E; % Easting (in meters)
y0 = N; % Northing (in meters)
z0 = Z; % Ship altitude

ship_trajectory = [x0 y0 z0];

GroundRange_port = zeros(size(port_ang));
Depth_port       = zeros(size(port_ang));
GroundRange_stbd = zeros(size(stbd_ang));
Depth_stbd       = zeros(size(stbd_ang));
for i = 1:nbPing67
    iPing=idxPing67(i)
    
    fs=Ping(iPing).truepix.H0_RxSampleRate;
    c=Ping(iPing).truepix.H0_SoundSpeed;
    
    range=[1:double(Ping(iPing).truepix.D(1).TotalSamples)]*double(c/2/fs);
    
    iSample=1:double(Ping(iPing).truepix.D(1).TotalSamples);%find(port_mag(i,:)>mag_thres);
    GroundRange_port(i,:) = range(iSample).*sin(port_ang(i,iSample)-roll67(i)*pi/180);
    Depth_port(i,:)       = range(iSample).*cos(port_ang(i,iSample)-roll67(i)*pi/180);
    
    GroundRange_stbd(i,:) = range(iSample).*sin(stbd_ang(i,iSample)-roll67(i)*pi/180);
    Depth_stbd(i,:)       = range(iSample).*cos(stbd_ang(i,iSample)-roll67(i)*pi/180);
end

GroundRange=[GroundRange_port GroundRange_stbd];
Depth=[Depth_port Depth_stbd];

alpha  = heading67+yaw67; % Heading angle + Yaw

points = zeros(size(Depth,1), size(Depth,2), 4);
delta_x = zeros(size(Depth,1),size(Depth,2));
delta_y = zeros(size(Depth,1),size(Depth,2));

delta_d_x = zeros(size(Depth,1),size(Depth,2));
delta_d_y = zeros(size(Depth,1),size(Depth,2));
delta_d_z = zeros(size(Depth,1),size(Depth,2));
for iPing = 1:size(points,1)
    for iSample = 1:size(points,2)
        
        delta_x(iPing,iSample) = GroundRange(iPing,iSample)*sind(alpha(iPing)+90); % Angle of the swath (perpendicular to the Heading angle)
        delta_y(iPing,iSample) = GroundRange(iPing,iSample)*cosd(alpha(iPing)+90);
        
        %% Pitch correction
        delta_d_z(iPing,iSample) = 0; %Depth(iPing,iSample)*cosd(pitch67(iPing));
        delta_d_x(iPing,iSample) = 0; %Depth(iPing,iSample)*sind(pitch67(iPing))*cosd(alpha(iPing));
        delta_d_y(iPing,iSample) = 0; %Depth(iPing,iSample)*sind(pitch67(iPing))*sind(alpha(iPing));
        
        %% (X, Y, Z)
        points(iPing,iSample,1) = x0(iPing) + delta_x(iPing,iSample) + delta_d_x(iPing,iSample);   % X-coordinate
        points(iPing,iSample,2) = y0(iPing) + delta_y(iPing,iSample) + delta_d_y(iPing,iSample);   % Y-coordinate
        points(iPing,iSample,3) = z0(iPing) - Depth(iPing,iSample) - delta_d_z(iPing,iSample);     % Z-coordinate
    end
    %% Add Intensity (from Magnitude Re-quantified for visual purpose)
    points(iPing,:,4) = [port_mag_req(iPing,:) stbd_mag_req(iPing,:)];
end

%%
figure
pcshow(points(:,:,1:3))




%%
nb_pts = size(Depth,1)*size(Depth,2); % to transform matrix to 1-D array
points = reshape(points, [nb_pts,4]);

mag_thres=200;
points_new=[];
for i = 1:nbPing67
    iPing=idxPing67(i)
    iSample_port=find(port_mag(i,:)>mag_thres);
    iSample_stbd=find(stbd_mag(i,:)>mag_thres);
    
    pts_port=reshape(points(i,iSample_port,:),[length(iSample_port), 4]);
    pts_stbd=reshape(points(i,size(port_ang,2)+iSample_stbd,:),[length(iSample_stbd), 4]);
    points_new=[points_new; pts_port; pts_stbd];
end


figure(idxFile)
pcshow(points_new(:,1:3),points_new(:,4))
hold on
plot3(E,N,Z,'r>'), grid on
legend('3D points','Sonar position')
xlabel('Easting (UTM)')
ylabel('Northing (UTM)')
zlabel('Altitude')


%%

mag_thres=50;
points_thres_700=[];
for i = idx700
    iPing=idxPing67(i)
    iSample_port=find(port_mag(i,:)>mag_thres);
    iSample_stbd=find(stbd_mag(i,:)>mag_thres);
    
    pts_port=reshape(points(i,iSample_port,:),[length(iSample_port), 4]);
    pts_stbd=reshape(points(i,size(port_ang,2)+iSample_stbd,:),[length(iSample_stbd), 4]);
    points_thres_700=[points_thres_700; pts_port; pts_stbd];
end

points_thres_170=[];
for i = idx170
    iPing=idxPing67(i)
    iSample_port=find(port_mag(i,:)>mag_thres);
    iSample_stbd=find(stbd_mag(i,:)>mag_thres);
    
    pts_port=reshape(points(i,iSample_port,:),[length(iSample_port), 4]);
    pts_stbd=reshape(points(i,size(port_ang,2)+iSample_stbd,:),[length(iSample_stbd), 4]);
    points_thres_170=[points_thres_170; pts_port; pts_stbd];
end

points_thres_330=[];
for i = idx330
    iPing=idxPing67(i)
    iSample_port=find(port_mag(i,:)>mag_thres);
    iSample_stbd=find(stbd_mag(i,:)>mag_thres);
    
    pts_port=reshape(points(i,iSample_port,:),[length(iSample_port), 4]);
    pts_stbd=reshape(points(i,size(port_ang,2)+iSample_stbd,:),[length(iSample_stbd), 4]);
    points_thres_330=[points_thres_330; pts_port; pts_stbd];
end


points_thres_450=[];
for i = idx450
    iPing=idxPing67(i)
    iSample_port=find(port_mag(i,:)>mag_thres);
    iSample_stbd=find(stbd_mag(i,:)>mag_thres);
    
    pts_port=reshape(points(i,iSample_port,:),[length(iSample_port), 4]);
    pts_stbd=reshape(points(i,size(port_ang,2)+iSample_stbd,:),[length(iSample_stbd), 4]);
    points_thres_450=[points_thres_450; pts_port; pts_stbd];
end

points170=points(idx170,:,:);
points330=points(idx330,:,:);
points450=points(idx450,:,:);
points700=points(idx700,:,:);

figure
pcshow(points330(:,:,1:3))

figure
pcshow(points_thres_330(:,1:3),points_thres_330(:,4))

figure
pcshow(points_thres_700(:,1:3),points_thres_700(:,4))
hold on
pcshow(points_thres_330(:,1:3),'r')

figure
subplot(221), pcshow(points_thres_170(:,1:3),points_thres_170(:,4))
xlabel('Port side Magnitude (170 kHz)')
subplot(222), pcshow(points_thres_330(:,1:3),points_thres_330(:,4)),
xlabel('Port side Magnitude (330 kHz)')
subplot(223), pcshow(points_thres_450(:,1:3),points_thres_450(:,4)), 
xlabel('Port side Magnitude (450 kHz)')
subplot(224), pcshow(points_thres_700(:,1:3),points_thres_700(:,4)), 
xlabel('Port side Magnitude (700 kHz)')


%%
points_normal_new=[];
for i = 1:length(idx_normal)
    iPing=idxPing67(idx_normal(i))
    iSample_port=find(port_mag(idx_normal(i),:)>mag_thres);
    iSample_stbd=find(stbd_mag(idx_normal(i),:)>mag_thres);
    
    pts_port=reshape(points_normal(i,iSample_port,:),[length(iSample_port), 4]);
    pts_stbd=reshape(points_normal(i,size(port_ang,2)+iSample_stbd,:),[length(iSample_stbd), 4]);
    points_normal_new=[points_normal_new; pts_port; pts_stbd];
end


figure(idxFile)
pcshow(points_normal_new(:,1:3),points_normal_new(:,4))
hold on
plot3(E,N,Z,'r>'), grid on
legend('3D points','Sonar position')
xlabel('Easting (UTM)')
ylabel('Northing (UTM)')
zlabel('Altitude')

%%


points_wc_new=[];
for i = 1:length(idx)
    iPing=idxPing67(idx(i))
    iSample_port=find(port_mag(idx(i),:)>mag_thres);
    iSample_stbd=find(stbd_mag(idx(i),:)>mag_thres);
    
    pts_port=reshape(points_wc(i,iSample_port,:),[length(iSample_port), 4]);
    pts_stbd=reshape(points_wc(i,size(port_ang,2)+iSample_stbd,:),[length(iSample_stbd), 4]);
    points_wc_new=[points_wc_new; pts_port; pts_stbd];
end


figure(idxFile)
pcshow(points_wc_new(:,1:3),points_wc_new(:,4))
hold on
plot3(E,N,Z,'r>'), grid on
legend('3D points','Sonar position')
xlabel('Easting (UTM)')
ylabel('Northing (UTM)')
zlabel('Altitude')