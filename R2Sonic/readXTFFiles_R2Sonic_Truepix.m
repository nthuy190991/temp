% Created by TH Nguyen, ULaval
% Last modified on: 2019/11/07


%% Parameters
v=1; % verbose mode

currentFolder = pwd;

data_dir = '/Users/thanhhuynguyen/Storage/r2SonicChallenge19/XTF_truepix/'
dirResult = '/Users/thanhhuynguyen/Storage/r2SonicChallenge19/XTF_truepix/mat/'
% data_dir = '/Users/thanhhuynguyen/Storage/r2SonicChallenge19/Eastcom_20190529 100-200-400kHz_stationary/'
% dirResult = '/Users/thanhhuynguyen/Storage/r2SonicChallenge19/Eastcom_20190529 100-200-400kHz_stationary/'
% data_dir = '/Users/thanhhuynguyen/Storage/r2SonicChallenge19/DC_20190226/'

% data_dir = '/Users/thanhhuynguyen/Storage/DATA/02_DiamondCreek/DC_20190226_truepix/Run D/'
% dirResult = '/Users/thanhhuynguyen/Storage/DATA/02_DiamondCreek/mat_truepix/Run D/'

cd(data_dir)
list = dir('*.xtf');

cd(currentFolder)

%% Reading XTF
for idxFile = 1:size(list,1)
    
    clear xtffile_header ChanInfo Ping
    fileName = list(idxFile).name

    %% Program body

    fileXTF  = [data_dir, fileName];
    fileID   = fopen(fileXTF);
%     PingType = 'Klein';

    %% XTF File header.
    %  Total of 1024 bytes.

    xtffile_header.FileFormat              = fread(fileID,1);
    xtffile_header.SystemType              = fread(fileID,1);
    xtffile_header.RecordingProgramName    = fread(fileID,[1 8],'*char');
    xtffile_header.RecordingProgramVersion = fread(fileID,[1 8],'*char');
    xtffile_header.SonarName               = fread(fileID,[1 16],'*char');
    xtffile_header.SonarType               = fread(fileID,1,'uint16');  
    xtffile_header.NoteString              = fread(fileID,[1 64],'*char');
    xtffile_header.ThisFileName            = fread(fileID,[1 64],'*char');
    xtffile_header.NavUnits                = fread(fileID,1,'uint16');
    xtffile_header.NumberOfSonarChannels   = fread(fileID,1,'uint16');
    xtffile_header.NumberOfBathymetryChannels      = fread(fileID,1,'uint16');
    xtffile_header.NumberOfForwardLookArrays       = fread(fileID,1,'uint16');
    xtffile_header.NumberOfEchoStrengthChannels    = fread(fileID,1,'uint16');
    xtffile_header.NumberOfInterferometriyChannels = fread(fileID,1,'*char');
    xtffile_header.Reserved1               = fread(fileID,1,'*char');
    xtffile_header.Reserved2               = fread(fileID,1,'uint16');
    xtffile_header.ReferencePointHeight    = fread(fileID,1,'float');

    %% Navigation System Parameters
%     xtffile_header.SystemSerialNumber      = fread(fileID,2,'*char');
    xtffile_header.ProjectionType          = fread(fileID,[1 12],'*char');
    xtffile_header.SpheriodType            = fread(fileID,[1 10],'*char');
    xtffile_header.NavigationLatency       = fread(fileID,1,'long');
    xtffile_header.OriginX                 = fread(fileID,1,'float');
    xtffile_header.OriginY                 = fread(fileID,1,'float');
    xtffile_header.NavOffsetY              = fread(fileID,1,'float');
    xtffile_header.NavOffsetX              = fread(fileID,1,'float');
    xtffile_header.NavOffsetZ              = fread(fileID,1,'float');
    xtffile_header.NavOffsetYaw            = fread(fileID,1,'float');
    xtffile_header.MRUOffsetY              = fread(fileID,1,'float');
    xtffile_header.MRUOffsetX              = fread(fileID,1,'float');
    xtffile_header.MRUOffsetZ              = fread(fileID,1,'float');
    xtffile_header.MRUOffsetYaw            = fread(fileID,1,'float');
    xtffile_header.MRUOffsetPitch          = fread(fileID,1,'float');
    xtffile_header.MRUOffsetRoll           = fread(fileID,1,'float');

    %% Channel information CHANINFO

    for iChannel=1:6
        xtffile_header.ChanInfo(iChannel).TypeOfChannel    = fread(fileID,1);
        xtffile_header.ChanInfo(iChannel).SubChannelNumber = fread(fileID,1);
        xtffile_header.ChanInfo(iChannel).CorrectionFlags  = fread(fileID,1,'uint16');
        xtffile_header.ChanInfo(iChannel).UniPolar         = fread(fileID,1,'uint16');
        xtffile_header.ChanInfo(iChannel).BytesPerSample   = fread(fileID,1,'uint16');
        xtffile_header.ChanInfo(iChannel).Reserved         = fread(fileID,1,'uint32');
        xtffile_header.ChanInfo(iChannel).ChannelName      = fread(fileID,[1 16],'*char');
        xtffile_header.ChanInfo(iChannel).VoltScale        = fread(fileID,1,'float');
        xtffile_header.ChanInfo(iChannel).Frequency        = fread(fileID,1,'float');
        xtffile_header.ChanInfo(iChannel).HorizBeamAngle   = fread(fileID,1,'float');
        xtffile_header.ChanInfo(iChannel).TiltAngle        = fread(fileID,1,'float');
        xtffile_header.ChanInfo(iChannel).BeamWidth        = fread(fileID,1,'float');
        xtffile_header.ChanInfo(iChannel).OffsetX          = fread(fileID,1,'float');
        xtffile_header.ChanInfo(iChannel).OffsetY          = fread(fileID,1,'float');
        xtffile_header.ChanInfo(iChannel).OffsetZ          = fread(fileID,1,'float');
        xtffile_header.ChanInfo(iChannel).OffsetYaw        = fread(fileID,1,'float');
        xtffile_header.ChanInfo(iChannel).OffsetPitch      = fread(fileID,1,'float');
        xtffile_header.ChanInfo(iChannel).OffsetRoll       = fread(fileID,1,'float');
        xtffile_header.ChanInfo(iChannel).BeamsPerArray    = fread(fileID,1,'uint16');
        xtffile_header.ChanInfo(iChannel).Latency          = fread(fileID,1,'float');
        xtffile_header.ChanInfo(iChannel).ReservedArea2    = fread(fileID,[1 50],'*char');
    end

    %% Ping Klein
    iPing = 0;
    while (1)
        iPing=iPing+1;

        Ping(iPing).MagicNumber = fread(fileID,1,'uint16');
        if (feof(fileID))
            break
        end
        if (Ping(iPing).MagicNumber ~= 64206)
            display(['Wrong Magic Number at Ping ' num2str(iPing)])
        end

        Ping(iPing).HeaderType = fread(fileID,1);

        display(['No Ping: ' num2str(iPing) ', HeaderType: ' num2str(Ping(iPing).HeaderType)])
        switch Ping(iPing).HeaderType
            case 0
                %% Read XTF Ping Header
                Ping(iPing).xtfping_header = XTFPINGHEADER(fileID);
                if (Ping(iPing).xtfping_header.NumChansToFollow>0)
                    %% Read XTF Ping Channel Header
                    Ping(iPing).xtfpingchan_header1 = XTFPINGCHANHEADER(fileID);
%                     Ping(iPing).chan1Sample = fread(fileID,[Ping(iPing).xtfpingchan_header1.NumSamples 1],'uint16');
                    Ping(iPing).chan1Sample = fread(fileID,[1024 1], 'uint16');
                end
                if (Ping(iPing).xtfping_header.NumChansToFollow>1)
                    Ping(iPing).xtfpingchan_header2 = XTFPINGCHANHEADER(fileID);
%                     Ping(iPing).chan2Sample = fread(fileID,[Ping(iPing).xtfpingchan_header1.NumSamples 1],'uint16');
                    Ping(iPing).chan2Sample = fread(fileID,[1024 1], 'uint16');
                end
%                 if (Ping(iPing).xtfping_header.NumChansToFollow>2)
%                     Ping(iPing).xtfpingchan_header3 = XTFPINGCHANHEADER(fileID);
%                     Ping(iPing).chan3Sample = fread(fileID,[Ping(iPing).xtfpingchan_header3.NumSamples 1],'uint16');
%                 end
%                 if (Ping(iPing).xtfping_header.NumChansToFollow>3)
%                     Ping(iPing).xtfpingchan_header4 = XTFPINGCHANHEADER(fileID);
%                     Ping(iPing).chan4Sample = fread(fileID,[Ping(iPing).xtfpingchan_header4.NumSamples 1],'uint16');
%                 end
% 
%                 skip=fread(fileID,48);  % padding pas elegant KLEIN
            case 2
                %% Read XTF Bathy data packet
                Ping(iPing).xtfbathy_header = XTFPINGHEADER(fileID);
                Ping(iPing).r_theta = R_THETA_DATA_OLD(fileID);
                nbByteToSkip = Ping(iPing).xtfbathy_header.NumBytesThisRecord-256-36-Ping(iPing).r_theta.beam_count*2-floor(Ping(iPing).r_theta.beam_count/2)-2;
                skip = fread(fileID, nbByteToSkip);
                
            case 3
                %% Read XTF Attitude data packet
                Ping(iPing).attitudedata = XTFATTITUDEDATA(fileID);
            case 6
                Ping(iPing).raw_serial = XTFRAWSERIALHEADER(fileID);
            case 8 
                %% Read XTF Hidden Ping Header
                Ping(iPing).hidden_ping_header = XTFPINGHEADER(fileID);
                if (Ping(iPing).hidden_ping_header.NumChansToFollow>0)
                    %% Read XTF Ping Channel Header
                    Ping(iPing).xtfpingchan_header1 = XTFPINGCHANHEADER(fileID);
                    Ping(iPing).chan1Sample = fread(fileID,[Ping(iPing).xtfpingchan_header1.NumSamples 1],'uint16');
                end
                if (Ping(iPing).hidden_ping_header.NumChansToFollow>1)
                    Ping(iPing).xtfpingchan_header2 = XTFPINGCHANHEADER(fileID);
                    Ping(iPing).chan2Sample = fread(fileID,[Ping(iPing).xtfpingchan_header1.NumSamples 1],'uint16');
                end
                if (Ping(iPing).hidden_ping_header.NumChansToFollow>2)
                    Ping(iPing).xtfpingchan_header3 = XTFPINGCHANHEADER(fileID);
                    Ping(iPing).chan3Sample = fread(fileID,[Ping(iPing).xtfpingchan_header3.NumSamples 1],'uint16');
                end
                if (Ping(iPing).hidden_ping_header.NumChansToFollow>3)
                    Ping(iPing).xtfpingchan_header4 = XTFPINGCHANHEADER(fileID);
                    Ping(iPing).chan4Sample = fread(fileID,[Ping(iPing).xtfpingchan_header4.NumSamples 1],'uint16');
                end
                skip = fread(fileID, 48);  % padding pas elegant KLEIN

            case 15
                Ping(iPing).xtf_header_highspeed_sensor2 = XTF_HEADER_HIGHSPEED_SENSOR(fileID);

            case 42
                Ping(iPing).isis_navigation = ISIS_NAVIGATION(fileID);
            
            case 67 
                %% R2SONIC TRUEPIX
                Ping(iPing).xtfbathy_header = XTFPINGHEADER(fileID);
                Ping(iPing).truepix = R2SONIC_TRUEPIX(fileID);
                
                sizeD = 344; % size of each UDP-splitted data section
                nbByteToSkip = Ping(iPing).xtfbathy_header.NumBytesThisRecord-256-Ping(iPing).truepix.PacketSize-sizeD*size(Ping(iPing).truepix.D,2);
                % (256: nb of bytes of the bathyheader)
                skip = fread(fileID, nbByteToSkip);
                
                
            case 73
                %% Read XTF EdgeTech data packet
                Ping(iPing).xtfping_header = XTFPINGHEADER(fileID);
                Ping(iPing).xtfping_ET4600 = XTFPINGBATHY_ET4600(fileID);

            case 84
                Ping(iPing).isis_gyro = ISIS_GYRO(fileID);
                
            case 107
                %% POSRAW_NAVIGATION
                Ping(iPing).xtf_posraw_navigation = XTFPOSRAWNAVIGATION(fileID);
            case 200 
                %% XTF_HEADER_USERDEFINED
                Ping(iPing).xtf_struct_userdefined = XTF_STRUCT_USER_DEFINED(fileID);

            otherwise
                display('Not recognized HeaderType')
        end

    end
    nbPing = size(Ping,2)-1;
    Ping = Ping(1:nbPing);
    fclose(fileID);

    save(strcat(dirResult,fileName,'.mat'),'Ping','xtffile_header')
end

