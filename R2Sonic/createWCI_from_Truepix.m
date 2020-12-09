%% Load xtf mat file
close all
clc

verbose=0;
currentFolder = pwd;
data_dir = '/Users/thanhhuynguyen/Storage/r2SonicChallenge19/XTF_truepix/'
dirResult = '/Users/thanhhuynguyen/Storage/r2SonicChallenge19/XTF_truepix/'
% data_dir = '/Users/thanhhuynguyen/Storage/r2SonicChallenge19/Eastcom_20190529 100-200-400kHz_stationary/'
% dirResult = '/Users/thanhhuynguyen/Storage/r2SonicChallenge19/Eastcom_20190529 100-200-400kHz_stationary/'

cd(data_dir)
list = dir('*.xtf');

cd(currentFolder)
idxFile=1
fileName = list(idxFile).name
load(strcat(dirResult,'mat/',fileName,'.mat'))

iPing=41;
d_all_packet=Ping(iPing).truepix.D_all_packets;
angle_factor=double(Ping(iPing).truepix.D(1).AngleScalingFactor);
C=10;
data=[C*log(1+d_all_packet(:,1)*2^8/2^16) d_all_packet(:,2)*angle_factor;
      C*log(1+d_all_packet(:,3)*2^8/2^16) d_all_packet(:,4)*angle_factor];

figure
plot(data(:,1))


%   figure
%     subplot(221), plot(Ping(45).truepix.D_all_packets(:,1)), grid on, xlabel('Samples'), ylabel('Port Mag. [\muPa]')
%     subplot(223), plot(Ping(45).truepix.D_all_packets(:,2).*Ping(45).truepix.D(1).AngleScalingFactor), grid on, xlabel('Samples'), ylabel('Port Angle [rad]')
%     subplot(222), plot(Ping(45).truepix.D_all_packets(:,3)), grid on, xlabel('Samples'), ylabel('Stbd Mag. [\muPa]')
%     subplot(224), plot(Ping(45).truepix.D_all_packets(:,4).*Ping(45).truepix.D(1).AngleScalingFactor), grid on, xlabel('Samples'), ylabel('Stbd Angle [rad]')
 

x=linspace(-20,20,500); %across distance
y=linspace(0,12,500); %depth

wci=zeros(length(y),length(x));

fs=Ping(iPing).truepix.H0_RxSampleRate;
c=Ping(iPing).truepix.H0_SoundSpeed;
total_samples=double(Ping(iPing).truepix.D(1).TotalSamples);
range=[1:total_samples]*double(c/2/fs);
     
xi=zeros(1,length(data));
yi=zeros(1,length(data));
for i=1:length(data)
    if i>total_samples
        r=range(i-total_samples);
    else
        r=range(i);
    end
    theta=data(i,2);
    xi(i)=r*sin(theta);
    yi(i)=r*cos(theta);
    
    idxx=min(find(x>=xi(i)));
    idxy=min(find(y>=yi(i)));
    wci(idxy,idxx)=data(i,1);
end

figure
imagesc(x,y,wci)
xlabel('Across-track distance [m]')
ylabel('Depth [m]')
% axis xy

figure
scatter(xi,yi,20,data(:,1),'filled')
grid on
axis ij
% hold on
% for i=1:length(data)
%     plot(xi(i),yi(i),'.')
% end

%%

wci_acc=zeros(length(y),length(x));

for iPing=idxPing67
d_all_packet=Ping(iPing).truepix.D_all_packets;
angle_factor=double(Ping(iPing).truepix.D(1).AngleScalingFactor);
C=10;
data=[C*log(1+d_all_packet(:,1)*2^8/2^16) d_all_packet(:,2)*angle_factor;
      C*log(1+d_all_packet(:,3)*2^8/2^16) d_all_packet(:,4)*angle_factor];

fs=Ping(iPing).truepix.H0_RxSampleRate;
c=Ping(iPing).truepix.H0_SoundSpeed;
total_samples=double(Ping(iPing).truepix.D(1).TotalSamples);
range=[1:total_samples]*double(c/2/fs);
     
xi=zeros(1,length(data));
yi=zeros(1,length(data));
for i=1:length(data)
    if i>total_samples
        r=range(i-total_samples);
    else
        r=range(i);
    end
    theta=data(i,2);
    xi(i)=r*sin(theta);
    yi(i)=r*cos(theta);
    
    idxx=min(find(x>=xi(i)));
    idxy=min(find(y>=yi(i)));
    wci_acc(idxy,idxx)=wci_acc(idxy,idxx)+data(i,1);
end
end

figure
imagesc(x,y,wci_acc)
xlabel('Across-track distance [m]')
ylabel('Depth [m]')
% 
% grid on
% myColorMap = parula(256);
% myColorMap(1,:) = 1;
% colormap(myColorMap);
% colorbar