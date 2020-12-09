% Created by TH Nguyen, ULaval
% Last modified on: 2019/11/07

%% Read D0 or D1 Truepix data in one UDP packet
function D=R2SONIC_TRUEPIX_D(fileID, D0_or_D1)

D.SectionName = D0_or_D1;                                      % 'D0'
D.SectionSize = swapbytes(uint16(fread(fileID,1,'uint16')));   % [bytes] size ofthis entire section
D.PingNumber = swapbytes(uint32(fread(fileID,1,'uint32')));    % pings since power-up or reboot
D.TotalSamples = swapbytes(uint32(fread(fileID,1,'uint32')));  % number of samples in entire time series (sample rate is H0_RxSampleRate)
D.FirstSample = swapbytes(uint32(fread(fileID,1,'uint32')));   % first sample of this section relative to zero range

D.Samples = swapbytes(uint16(fread(fileID,1,'uint16')));       % number of samples in this section
D.reserved = swapbytes(uint16(fread(fileID,1,'uint16')));      % reserved for future use

D.MagnitudeScaling = swapbytes(single(fread(fileID,[1 8],'float32')));  % to be determined, 0=ignore

if strcmp(D0_or_D1,'D0')
    %% section D0: 16-bit magnitude data (present only during "magnitude only" mode)
    % disp('D0')
    D.Data=zeros(D.Samples,2);
    for i=1:D.Samples
        D.Data(i,1) = swapbytes(uint16(fread(fileID,1,'uint16'))); % PortMagnitude, [micropascals] = PortMagnitude * (tbd function of sample number and D0_MagnitudeScaling[8])
        D.Data(i,2) = swapbytes(uint16(fread(fileID,1,'uint16'))); % StbdMagnitude, similar but starboard side
    end
    
elseif strcmp(D0_or_D1,'D1')
    %% section D1: 16-bit magnitude and direction data (present only during "magnitude+direction" mode)
    D.AngleScalingFactor = swapbytes(single(fread(fileID,1,'float32')));
    D.Data=zeros(D.Samples,4);
    for i=1:D.Samples
        D.Data(i,1) = swapbytes(uint16(fread(fileID,1,'uint16'))); % PortMagnitude; // [micropascals] = PortMagnitude * (tbd function of sample number and D1_MagnitudeScaling[8])
        D.Data(i,2) = swapbytes(int16(fread(fileID,1,'int16')));   % [radians from array centerline (positive towards starboard)] = PortAngle * D1_AngleScalingFactor
        D.Data(i,3) = swapbytes(uint16(fread(fileID,1,'uint16'))); % similar but starboard side
        D.Data(i,4) = swapbytes(int16(fread(fileID,1,'int16')));   % similar but starboard side
    end
end