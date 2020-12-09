dir_data='C:\Users\thnguyen\Dropbox\Thesis\Projet_R2Sonic\samples\';
filename='CRabvDC-SUSP-SED_1997-2019.xlsx';
T = readtable(strcat(dir_data,filename),'FileType','spreadsheet')

%% Get fields 

DISCHARGE_COL=T.Q_cfs_;
DATE_COL=T.date;
TIME_COL=T.startTimeOr;
% To convert to hours: 0.3333*24 = 8 (hours)
% To convert to minutes: 0.3333*24 *60 = 480 (minutes)
% To convert to seconds: 0.3333*24 *60*60 = 28800 (seconds)

SAMPLING_METHOD_COL=T.SamplingMethod;

WATER_DEPTH_COL=T.Water;
SAMPLE_ELE_COL=T.Sample;

SILT_CLAY_CONC_LAB_COL=T.Laboratory_determinedValues;
SAND_CONC_LAB_COL=T.Var28;
SAND_D50_LAB_COL=T.Var30;


%%
SILT_CLAY_CONC_LAB=zeros(length(SILT_CLAY_CONC_LAB_COL)-5+1,1);
SAND_CONC_LAB=zeros(length(SAND_CONC_LAB_COL)-5+1,1);
SAND_D50_LAB=zeros(length(SAND_D50_LAB_COL)-5+1,1);

DATE={};
SAMPLING_METHOD={};
TIMESTAMP=zeros(length(SAND_D50_LAB_COL)-5+1,1);

k=1;
for i=5:length(SILT_CLAY_CONC_LAB_COL)
%     DATE(k)=strcat(DATE_COL(i),' ',TIME_COL(i));
    s=DATE_COL{i};
    s2=strsplit(s,'/');
    DATE(k)={strcat(s2{1},'/',s2{3}(3:4))}; % format MM/YY
    
    time_in_sec=str2double(TIME_COL{i})*24*60*60;
    t=datetime(str2double(s2{3}),str2double(s2{1}),str2double(s2{2}));
    TIMESTAMP(k)=posixtime(t)+time_in_sec;
    
    SAMPLING_METHOD(k)=SAMPLING_METHOD_COL(i);
    
    SILT_CLAY_CONC_LAB(k)=str2double(SILT_CLAY_CONC_LAB_COL(i));
    SAND_CONC_LAB(k)=str2double(SAND_CONC_LAB_COL(i));
    SAND_D50_LAB(k)=str2double(SAND_D50_LAB_COL(i));
    k=k+1;
end


% Let the first sample be the time origin
TIMESTAMP=TIMESTAMP-TIMESTAMP(1);


%% Concentration plots colored by sampling method
figure
hold on
color='rgbmcyk';
marker_SC='xxxxxxx';%'.+xosv^';
for i=1:length(SILT_CLAY_CONC_LAB)
    if strcmp(SAMPLING_METHOD{i},'EDI')
        p1=plot(i,SILT_CLAY_CONC_LAB(i),'Marker',marker_SC(1),'Color',color(1));
    elseif strcmp(SAMPLING_METHOD{i},'EWI')
        p2=plot(i,SILT_CLAY_CONC_LAB(i),'Marker',marker_SC(2),'Color',color(2));
    elseif strcmp(SAMPLING_METHOD{i},'dip')
        p3=plot(i,SILT_CLAY_CONC_LAB(i),'Marker',marker_SC(3),'Color',color(3));
    elseif strcmp(SAMPLING_METHOD{i},'individual vertical of EDI')
        p4=plot(i,SILT_CLAY_CONC_LAB(i),'Marker',marker_SC(4),'Color',color(4));
    elseif strcmp(SAMPLING_METHOD{i},'point')
        p5=plot(i,SILT_CLAY_CONC_LAB(i),'Marker',marker_SC(5),'Color',color(5));
    elseif strcmp(SAMPLING_METHOD{i},'pump')
        p6=plot(i,SILT_CLAY_CONC_LAB(i),'Marker',marker_SC(6),'Color',color(6));
    elseif strcmp(SAMPLING_METHOD{i},'single vertical')
        p7=plot(i,SILT_CLAY_CONC_LAB(i),'Marker',marker_SC(7),'Color',color(7));
    end
end
grid on
legend([p1,p2,p3,p4,p5,p6,p7],{'EDI','EWI','dip','individual vertical of EDI','point','pump','single vertical'})
set(gca,'xtick',[1:50:length(SILT_CLAY_CONC_LAB)],'xticklabel',DATE)


figure
hold on
color='rgbmcyk';
marker_S='.......';
for i=1:100%length(SILT_CLAY_CONC_LAB)
    if strcmp(SAMPLING_METHOD{i},'EDI')
        p1=plot(i,SAND_CONC_LAB(i),'Marker',marker_S(1),'Color',color(1));
    elseif strcmp(SAMPLING_METHOD{i},'EWI')
        p2=plot(i,SAND_CONC_LAB(i),'Marker',marker_S(2),'Color',color(2));
    elseif strcmp(SAMPLING_METHOD{i},'dip')
        p3=plot(i,SAND_CONC_LAB(i),'Marker',marker_S(3),'Color',color(3));
    elseif strcmp(SAMPLING_METHOD{i},'individual vertical of EDI')
        p4=plot(i,SAND_CONC_LAB(i),'Marker',marker_S(4),'Color',color(4));
    elseif strcmp(SAMPLING_METHOD{i},'point')
        p5=plot(i,SAND_CONC_LAB(i),'Marker',marker_S(5),'Color',color(5));
    elseif strcmp(SAMPLING_METHOD{i},'pump')
        p6=plot(i,SAND_CONC_LAB(i),'Marker',marker_S(6),'Color',color(6));
    elseif strcmp(SAMPLING_METHOD{i},'single vertical')
        p7=plot(i,SAND_CONC_LAB(i),'Marker',marker_S(7),'Color',color(7));
    end
end
grid on
legend([p1,p2,p3,p4,p5,p6,p7],{'EDI','EWI','dip','individual vertical of EDI','point','pump','single vertical'})
set(gca,'xtick',[1:50:length(SILT_CLAY_CONC_LAB)],'xticklabel',DATE)

%% Samples since 11/13/2013 (row=3177)

SILT_CLAY_CONC_LAB=zeros(length(SILT_CLAY_CONC_LAB_COL)-3177+1,1);
SAND_CONC_LAB=zeros(length(SAND_CONC_LAB_COL)-3177+1,1);
SAND_D50_LAB=zeros(length(SAND_D50_LAB_COL)-3177+1,1);

DISCHARGE=zeros(length(DISCHARGE_COL)-3177+1,1);

WATER_DEPTH=zeros(length(WATER_DEPTH_COL)-3177+1,1);
SAMPLE_ELE=zeros(length(SAMPLE_ELE_COL)-3177+1,1);

DATE={};
SAMPLING_METHOD={};
TIMESTAMP=zeros(length(SAND_D50_LAB_COL)-3177+1,1);

k=1;
for i=3177:length(SILT_CLAY_CONC_LAB_COL)
    i
    DISCHARGE(k)=str2double(DISCHARGE_COL(i));
    
%     DATE(k)=strcat(DATE_COL(i),' ',TIME_COL(i));
    s=DATE_COL{i};
    s2=strsplit(s,'/');
    DATE(k)={strcat(s2{1},'/',s2{3}(3:4))}; % format MM/YY
    
    if ~isnan(str2double(TIME_COL{i}))
        time_in_sec=str2double(TIME_COL{i})*24*60*60;
        t=datetime(str2double(s2{3}),str2double(s2{1}),str2double(s2{2}));
        TIMESTAMP(k)=posixtime(t)+time_in_sec;
    else
        s3=strsplit(TIME_COL{i},':');
        t=datetime(str2double(s2{3}),str2double(s2{1}),str2double(s2{2}),str2double(s3{1}),str2double(s3{2}),0);
        TIMESTAMP(k)=posixtime(t);
    end
    
    SAMPLING_METHOD(k)=SAMPLING_METHOD_COL(i);
    
    WATER_DEPTH(k)=str2double(WATER_DEPTH_COL(i));
    SAMPLE_ELE(k)=str2double(SAMPLE_ELE_COL(i));
    SILT_CLAY_CONC_LAB(k)=str2double(SILT_CLAY_CONC_LAB_COL(i));
    SAND_CONC_LAB(k)=str2double(SAND_CONC_LAB_COL(i));
    SAND_D50_LAB(k)=str2double(SAND_D50_LAB_COL(i));
    k=k+1;
end


% Let the first sample be the time origin
% TIMESTAMP=TIMESTAMP-TIMESTAMP(1);

figure
hold on
color='rgbmcyk';
marker_SC='xxxxxxx';%'.+xosv^';
for i=1:length(SILT_CLAY_CONC_LAB)
    if strcmp(SAMPLING_METHOD{i},'EDI')
        p1=plot(i,SILT_CLAY_CONC_LAB(i),'Marker',marker_SC(1),'Color',color(1));
    elseif strcmp(SAMPLING_METHOD{i},'EWI')
        p2=plot(i,SILT_CLAY_CONC_LAB(i),'Marker',marker_SC(2),'Color',color(2));
    elseif strcmp(SAMPLING_METHOD{i},'dip')
        p3=plot(i,SILT_CLAY_CONC_LAB(i),'Marker',marker_SC(3),'Color',color(3));
    elseif strcmp(SAMPLING_METHOD{i},'individual vertical of EDI')
        p4=plot(i,SILT_CLAY_CONC_LAB(i),'Marker',marker_SC(4),'Color',color(4));
    elseif strcmp(SAMPLING_METHOD{i},'point')
        p5=plot(i,SILT_CLAY_CONC_LAB(i),'Marker',marker_SC(5),'Color',color(5));
    elseif strcmp(SAMPLING_METHOD{i},'pump')
        p6=plot(i,SILT_CLAY_CONC_LAB(i),'Marker',marker_SC(6),'Color',color(6));
    elseif strcmp(SAMPLING_METHOD{i},'single vertical')
        p7=plot(i,SILT_CLAY_CONC_LAB(i),'Marker',marker_SC(7),'Color',color(7));
    end
end
grid on
legend([p1,p3,p4,p5,p7],{'EDI','dip','individual vertical of EDI','point','single vertical'})
% set(gca,'xtick',[1:length(SILT_CLAY_CONC_LAB)],'xticklabel',DATE)

figure
subplot(131), plot(WATER_DEPTH,'.'), grid on
subplot(132), plot(SAMPLE_ELE,'.'), grid on
subplot(133), plot(WATER_DEPTH+SAMPLE_ELE,'.'), grid on
