% Created by TH Nguyen, ULaval
% Last modified on: 2019/11/07

%% Read XTF ping header, reference: page 35 XTF_Format_X41

function xtf_posraw_navigation = XTFPOSRAWNAVIGATION(fileID)

% Read variables
xtf_posraw_navigation.SubChannelNumber = fread(fileID,1);
xtf_posraw_navigation.NumChansToFollow = fread(fileID,1,'uint16');
xtf_posraw_navigation.Reserved1 = fread(fileID,2,'uint16');
xtf_posraw_navigation.NumBytesThisRecord = fread(fileID,1,'uint32');
xtf_posraw_navigation.Year = fread(fileID,1,'uint16');
xtf_posraw_navigation.Month = fread(fileID,1);
xtf_posraw_navigation.Day = fread(fileID,1);
xtf_posraw_navigation.Hour = fread(fileID,1);
xtf_posraw_navigation.Minute = fread(fileID,1);
xtf_posraw_navigation.Second = fread(fileID,1);
xtf_posraw_navigation.MicroSeconds = fread(fileID,1,'uint16');

xtf_posraw_navigation.RawYcoordinate = fread(fileID,1,'double');
xtf_posraw_navigation.RawXcoordinate = fread(fileID,1,'double');
xtf_posraw_navigation.RawAltitude = fread(fileID,1,'double');

xtf_posraw_navigation.Pitch = fread(fileID,1,'float');
xtf_posraw_navigation.Roll = fread(fileID,1,'float');
xtf_posraw_navigation.Heave = fread(fileID,1,'float');
xtf_posraw_navigation.Heading = fread(fileID,1,'float');

xtf_posraw_navigation.Reserved2 = fread(fileID,1); 
end