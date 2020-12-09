% Created by TH Nguyen, ULaval
% Last modified on: 2019/11/14
%
% Calculation of backscattering strength per unit volume of ensonified
% water (Sv, in decibels)
% 
% References: 
% [1] Simmons et al. (2010). Monitoring Suspended Sediment Dynamics Using MBES. 
%     Journal of Hydraulic Engineering, 136(1), 45?49.
% [2] Gallaudet, T. C., & de Moustier, C. P. (2002). Multibeam volume acoustic 
%     backscatter imagery and reverberation measurements in the northeastern Gulf
%     of Mexico. The Journal of the Acoustical Society of America, 112(2), 489?503. 
%
% Change the index array in line 68 to display pings from different
% frequency, e.g. idx170 for 170KHz frequency pings.
%
%% Load .mat file
close all
clc

verbose=0;
currentFolder = pwd;

data_dir = '/02_DiamondCreek/mat_truepix/Run A/'
dirResult = data_dir

cd(data_dir)
list = dir('*.xtf.mat');


cd(currentFolder)

fileName = list(1).name
load(strcat(data_dir,fileName))

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


%% Display Sv as a function of pings of a chosen frequency in fan-shaped view

for i=idx450%idx700
    iPing=idxPing67(i);
    
    %% Get variables
    freq=Ping(iPing).truepix.H0_Frequency;
    SL0=Ping(iPing).truepix.H0_TxPower;
    
    pulse_length=Ping(iPing).truepix.H0_TxPulseWidth;
    c=Ping(iPing).truepix.H0_SoundSpeed;
    fs=Ping(iPing).truepix.H0_RxSampleRate;
    
    bw_vert=Ping(iPing).truepix.H0_TxBeamwidthVert;
    bw_horiz=Ping(iPing).truepix.H0_TxBeamwidthHoriz;
    
    steering_vert=Ping(iPing).truepix.H0_TxSteeringVert;
    steering_horiz=Ping(iPing).truepix.H0_TxSteeringHoriz;
    
    S_G=Ping(iPing).truepix.H0_RxSpreading;
    alpha_G=Ping(iPing).truepix.H0_RxAbsorption;
    G=Ping(iPing).truepix.H0_RxGain;
    
    if freq>5e5 % for 700 kHz
        swath=70/180*pi;
        alpha=175/1000;
        max_depth=15;
    else % for other frequencies
        swath=120/180*pi;
        if freq==17e4
            alpha=11/1000;
        elseif freq==33e4
            alpha=40/1000;
        else
            alpha=73/1000;
        end
        max_depth=20;
    end
    
    theta=-swath/2:bw_horiz*2:swath/2;
    delta_r=double(c*pulse_length/2);
    N=floor(max_depth/delta_r);
    r=(1:N)*delta_r;
    
    if ~exist('Sv_total','var')
        Sv_total=zeros(N-1,length(theta));
    end
    
    V=zeros(N,length(theta));
    for ir=1:N-1
        for itheta=1:length(theta)
            V(ir,itheta)=2/3*bw_vert*sin(bw_horiz/2)*(r(ir+1)^3-r(ir)^3);
        end
    end
    
    if (verbose)
        figure
        imagesc(V)
    end
    
    %% Begin calculation
    
    D=Ping(iPing).truepix.D_all_packets;
    port_mag=D(:,1)';
    stbd_mag=D(:,3)';
    all_mag=[fliplr(port_mag) stbd_mag];
    
    angle_factor=Ping(iPing).truepix.D(1).AngleScalingFactor;
    port_ang=D(:,2)'*angle_factor;
    stbd_ang=D(:,4)'*angle_factor;
    all_angle=[fliplr(port_ang) stbd_ang];
    
    C=20;
    port_mag_req=C*log10(1+port_mag(:,:)*2^8/2^16);
    stbd_mag_req=C*log10(1+stbd_mag(:,:)*2^8/2^16);
    all_mag_req=[fliplr(port_mag_req) stbd_mag_req];
    
    if (verbose)
        figure
        subplot(221), plot(fliplr(port_mag_req)), xlabel('Port side Magnitude')
        subplot(222), plot(stbd_mag_req), xlabel('Starboard side Magnitude')
        subplot(223), plot(fliplr(port_ang)), xlabel('Port side Angle')
        subplot(224), plot(stbd_ang), xlabel('Starboard side Angle')
    end
    
    total_samples=double(Ping(iPing).truepix.D(1).TotalSamples);
    range=[total_samples:-1:1 1:total_samples]*double(c/2/fs);

    RL=zeros(N,length(theta));
    for iSample=1:length(D)*2
        idx_r=min(find(r>range(iSample)));
        idx_theta=min(find(theta>all_angle(iSample)));
        
        RL(idx_r,idx_theta)=RL(idx_r,idx_theta)+all_mag_req(iSample);
    end
    if (verbose)
        figure
        imagesc(RL)
    end
    
    SL=zeros(N-1,length(theta));
    idx=min(find(r>1));
    SL(idx:end,:)=SL0; % TxPower [dB re 1 uPa] at 1 meter
    
    Sv=zeros(N-1,length(theta));
    for ir=idx:N-1
        R=r(ir);
        for itheta=1:length(theta)
            Sv(ir,itheta)=RL(ir,itheta) - SL(ir,itheta) + ...
                40*log10(R) + 2*alpha*R - 10*log10(V(ir,itheta))...
                -(S_G*log10(R) + 2*alpha_G*R/1000 + G);
        end
    end
    Sv_total=Sv_total+Sv;
    
    if (verbose)
        figure
        imagesc(Sv,[0 1000])
    end
    
    [rho, the] = meshgrid(-r(idx:end),theta);
    [X,Y] = pol2cart(the,rho);
    
    %% Display Sv
    
    figure(10)
    set(gcf, 'Position',  [100, 200, 500, 400])
    S = surf(Y,X,ones(size(X)));
    set(S,'FaceColor','Texturemap','CData',Sv(idx:end,:)','edgecolor','none')
    view(2)
    axis equal
    
    if freq>5e5
        title(['Frequency: 700000' ' - Ping ' num2str(iPing)])
    else
        title(['Frequency: ' num2str(freq) ' - Ping ' num2str(iPing)])
    end
    xlabel('Across distance [m]'), ylabel('Depth [m]')
    xlim([max_depth*sin(min(theta)) max_depth*sin(max(theta))])
    ylim([-max_depth 0])
%     caxis([-155 -130])
    caxis([-190 -150])
    colorbar
    
    %% Display nSv
    
%     figure(11)
%     set(gcf, 'Position',  [700, 200, 500, 400])
%     S = surf(Y,X,ones(size(X)));
%     set(S,'FaceColor','Texturemap','CData',mat2gray(Sv(idx:end,:))','edgecolor','none')
%     view(2)
%     axis equal
%     
%     if freq>5e5
%         title(['Frequency: 700000' ' - Ping ' num2str(iPing)])
%     else
%         title(['Frequency: ' num2str(freq) ' - Ping ' num2str(iPing)])
%     end
%     xlabel('Across distance [m]'), ylabel('Depth [m]')
%     xlim([max_depth*sin(min(theta)) max_depth*sin(max(theta))])
%     ylim([-max_depth 0])
%     
%     caxis([0 1])
%     colorbar
    
end


%% Display the all-ping-combined fan-shaped view of Sv


figure
hist(Sv_total(:),255)
disp('Use the histogram plot to determine the display range for the all-ping-combined view')

figure(12)
S = surf(Y,X,ones(size(X)));
set(S,'FaceColor','Texturemap','CData',Sv_total','edgecolor','none')
view(2)
axis equal

if freq>5e5
    title(['Frequency: 700000' ' - Ping ' num2str(iPing)])
else
    title(['Frequency: ' num2str(freq) ' - Ping ' num2str(iPing)])
end
xlabel('Across distance [m]'), ylabel('Depth [m]')
xlim([max_depth*sin(min(theta)) max_depth*sin(max(theta))])
ylim([-max_depth 0])

% caxis([-15.6e3 -15e3])
caxis([-4.65e4 -4.45e4])
colorbar

%% GWC

idx=4:4:length(tb_posX);
idx_normal=setdiff(1:length(tb_posX),idx);
for k=idx
    posX=tb_posX{k,:};
    posY=tb_posY{k,:};
    sampData=tb_sampData{k,:};
    BottomPick=tb_BottomPick{k,:};
    if plott
        figure(20)
        subplot(121)
        imagesc(sampData)
        
    subplot(122)
        [rho, the] = meshgrid(1:size(sampData,1),1:size(sampData,2));
%         [X,Y] = pol2cart(the,rho);

        S = surf(the,rho,ones(size(rho)));
        set(S,'FaceColor','Texturemap','CData',sampData','edgecolor','none')
        axis equal
    
    view(2)
        hold on
        
%         plot(BottomPick,'k.')
        
        hold off
        
%         caxis([-64 64])
        
%         c = colormap(jet);
        
%         c(1,:) = 1;
        
%         colormap(c)
        %             xlim([1 256])
        %             ylim([1 612])
        title(['Reading Ping Number: ' num2str(pidx)])
        drawnow
        
    end
    pidx=pidx+1;
    disp(['Reading Ping Number: ' num2str(pidx)])
end
