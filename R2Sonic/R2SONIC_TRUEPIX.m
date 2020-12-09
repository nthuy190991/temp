% Created by TH Nguyen, ULaval
% Last modified on: 2019/11/07

%% Read TruePix Data Format

% *** BEGIN TRUEPIX DATA FORMAT 0 ***
% TruePix is like sidescan with 3D relief. Each sonar ping produces a port
% and starboard time-series of data samples at the sonar's sample rate. Each
% sample contains the signal's magnitude (like sidescan) and across-track
% target direction angle (like bathymetry). After collecting many pings of
% data along a survey line, you now have a large array of data points with
% range, direction, and brightness. Apply noise reduction, and render the
% data as a textured 3D surface.
%  
% Two data formats are available: D0 provides magnitudes only, D1 provides
% magnitudes and direction angles. The GUI allows the user to choose the
% desired format.
%  
% The sonar generates one TruePix data set per ping. Each data set is
% usually split into multiple UDP packets. The D0 or D1 header includes
% FirstSample and Samples values to help you reassemble the full data set.
%  
% Someday you may be able to convert the 16-bit magnitude values to
% micropascals by applying a to-be-determined function involving the sample
% number and the MagnitudeScaling[] coefficients, but this conversion is not
% yet supported so these coefficients are zero. You can convert the direction
% angles from 16-bit values to radians by multiplying by AngleScalingFactor.


function truepix = R2SONIC_TRUEPIX(fileID)

% Note: This data packet is written in Big Endian. As Matlab uses Little 
% Endian by default, we need to convert it into Big Endian.
% In this case, fields in BYTE will be unchanged, but others will be added
% a swapbyte function in order to swap to Little Endian format

%% Read variables

truepix.PacketName = fread(fileID,[1 4],'*char');                     % 'TPX0'
truepix.PacketSize = swapbytes(uint32(fread(fileID,1,'uint32')));     % may be zero in UDP, otherwise:[bytes] size of this entire packet
truepix.DataStreamID = swapbytes(uint32(fread(fileID,1,'uint32')));   % reserved for future use

% Section H0: header (present only in first packet of each ping)
truepix.H0_SectionName = fread(fileID,[1 2],'*char');                 % 'H0'
truepix.H0_SectionSize = swapbytes(uint16(fread(fileID,1,'uint16'))); % [bytes] size ofthis entire section
truepix.H0_ModelNumber = fread(fileID,12,'uint8');                    % example "2024", unused chars are nulls
truepix.H0_SerialNumber = fread(fileID,12,'uint8');                   % example "100017", unused chars are nulls

truepix.H0_TimeSeconds = swapbytes(uint32(fread(fileID,1,'uint32')));     % [seconds] ping time relative to 0000 hours 1-Jan-1970, integer part
truepix.H0_TimeNanoseconds = swapbytes(uint32(fread(fileID,1,'uint32'))); % [nanoseconds] ping time relative to 0000 hours 1-Jan-1970, fraction part
truepix.H0_PingNumber = swapbytes(uint32(fread(fileID,1,'uint32')));      % pings since power-up or reboot

truepix.H0_PingPeriod = swapbytes(single(fread(fileID,1,'float32')));     % [seconds] time between most recent two pings
truepix.H0_SoundSpeed = swapbytes(single(fread(fileID,1,'float32')));     % [meters per second]
truepix.H0_Frequency = swapbytes(single(fread(fileID,1,'float32')));      % [hertz] sonar center frequency
truepix.H0_TxPower = swapbytes(single(fread(fileID,1,'float32')));        % [dB re 1 uPa at 1 meter]
truepix.H0_TxPulseWidth = swapbytes(single(fread(fileID,1,'float32')));   % [seconds]
truepix.H0_TxBeamwidthVert = swapbytes(single(fread(fileID,1,'float32')));   % [radians]
truepix.H0_TxBeamwidthHoriz = swapbytes(single(fread(fileID,1,'float32')));  % [radians]
truepix.H0_TxSteeringVert = swapbytes(single(fread(fileID,1,'float32')));    % [radians]
truepix.H0_TxSteeringHoriz = swapbytes(single(fread(fileID,1,'float32')));   % [radians]

truepix.H0_2026ProjTemp = swapbytes(uint16(fread(fileID,1,'uint16')));  % [hundredths of a degree Kelvin] 2026 projector temperature (divide value by 100, subtract 273.15 to get °C)

truepix.H0_VTX_Offset = swapbytes(int16(fread(fileID,1,'int16')));      % [hundredths of a dB] transmit voltage offset at time of ping (divide value by 100 to get dB)

truepix.H0_RxBandwidth = swapbytes(single(fread(fileID,1,'float32')));  % [hertz] 
truepix.H0_RxSampleRate = swapbytes(single(fread(fileID,1,'float32'))); % [hertz] sample rate of data acquisition and signal processing
truepix.H0_RxRange = swapbytes(single(fread(fileID,1,'float32')));      % user setting [meters]
truepix.H0_RxGain = swapbytes(single(fread(fileID,1,'float32')));       % user setting [multiply by 2 for dB]
truepix.H0_RxSpreading = swapbytes(single(fread(fileID,1,'float32')));  % [dB (times log range in meters)]
truepix.H0_RxAbsorption = swapbytes(single(fread(fileID,1,'float32'))); % [dB per kilometer]
truepix.H0_RxMountTilt = swapbytes(single(fread(fileID,1,'float32')));  % [radians]

truepix.H0_RxMiscInfo = swapbytes(uint32(fread(fileID,1,'uint32')));    % reserved for future use
truepix.H0_reserved = swapbytes(uint32(fread(fileID,1,'uint32')));      % reserved for future use

truepix.A0_MoreInfo_0 = swapbytes(single(fread(fileID,1,'float32')));   % 0 (reserved for future use)
truepix.A0_MoreInfo_1 = swapbytes(single(fread(fileID,1,'float32')));   % Z-offset, proj [metres]
truepix.A0_MoreInfo_2 = swapbytes(single(fread(fileID,1,'float32')));   % Y-offset, proj [metres]
truepix.A0_MoreInfo_3 = swapbytes(single(fread(fileID,1,'float32')));   % X-offset, proj [metres]
truepix.A0_MoreInfo_4 = swapbytes(single(fread(fileID,1,'float32')));   % 0 (reserved for future use)
truepix.A0_MoreInfo_5 = swapbytes(single(fread(fileID,1,'float32')));   % 0 (reserved for future use)
    
D0_or_D1=fread(fileID,[1 2],'*char');
truepix.D0_or_D1=strcmp(D0_or_D1,'D1');
truepix.D(1)=R2SONIC_TRUEPIX_D(fileID, D0_or_D1);

nb_samples=truepix.D(1).Samples;
D_all_packets=[truepix.D(1).Data]; % matrix that combines all packets

%% Combining UDP packets into one matrix (D_all_packets)
k=1;
while nb_samples<truepix.D(1).TotalSamples

    k=k+1;
    D0_or_D1=fread(fileID,[1 2],'*char');

    truepix.D0_or_D1=strcmp(D0_or_D1,'D1');
    truepix.D(k)=R2SONIC_TRUEPIX_D(fileID, D0_or_D1);
    
    D_all_packets=[D_all_packets; truepix.D(k).Data]; 

    nb_samples=nb_samples+truepix.D(k).Samples;
end

truepix.D_all_packets=D_all_packets;

% *** END TRUEPIX DATA FORMAT 0 ***

end