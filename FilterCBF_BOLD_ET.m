
% Written by Chiarelli Atonio Maria, University G. D'Annunzio of Chieti-Pescara
% Original Date: 03/15/2024
% Last UpDate: 03/15/2024



%The function filters ASL, BOLD and O2 and CO2 end tidals using butterworth
%digital filters


%Inputs:
%ASL structure, output of TE1DEXI2CBF.m
%BOLD structure, output of TE2DEXI2BOLD.m
% PETO2, pressure of end tidal O2 in mmHg, length of the signal equal to
% the length of BOLD.BOLDGM
% PETCO2, pressure of end tidal CO2 in mmHg, length of the signal equal to
% the length of BOLD.BOLDGM
%Parameters structure must report all required information 


%     Parameters.TR=4.4; %Time of echo of the BOLD and ASL signal, units of s
%     Parameters.THP=150; %high pass cut-off time, units of s
%     Parameters.TLP=20; %low pass cut-off time, units of s
%   Parameters.HPorder=4; % high-pass filter order
%   Parameters.LPorder=4; % low-pass filter order
%     Parameters.baseline=1:13; %baseline volumes
%     FigOn=1 creates matlab figures of the analysis outcome

%Outputs:
% fMRI_filt structure with the following field:
% CBF0, baseline CBF
% GM, GM binary mask as computed in TE1DEXI2CBF
% PaO20, baseline pressure of end tidal O2
% PaCO20, baseline pressure of end tidal CO2
% PETO2filt, end tidal O2 filtered
% PETCO2filt, end tidal CO2 filtered
% CBFfilt, filtered CBF signal
% BOLDfilt, filtered BOLD signal
% CBFGMfilt, median filtered CBF signal
% BOLDGMfilt, median filtered BOLD signal
% CBF0GM, median baseline CBF in the GM


function [fMRI_filt]=FilterCBF_BOLD_ET(ASL,BOLD,PETO2,PETCO2,Parameters,FigOn)

%% import the needed parameters
TR=Parameters.TR;
THP=Parameters.THP;
TLP=Parameters.TLP;
HPorder=Parameters.HPorder;
LPorder=Parameters.LPorder;
baseline=Parameters.baseline;



%% compute baseline values and changes in O2 and CO2
PaO20=mean(PETO2(baseline));
PaCO20=mean(PETCO2(baseline));
DCO2=PETCO2-PaCO20;
DO2=PETO2-PaO20;


%% filter end tidals
fsample=1/TR;
HP=1./THP;
LP=1./TLP;

if ~isempty(HP) && isempty(LP)
[NHP,DHP]=butter(HPorder,2*HP/fsample,'high');
DCO2f=filtfilt(NHP,DHP,DCO2);
DO2f=filtfilt(NHP,DHP,DO2);
end
if ~isempty(LP) && isempty(HP)
[NLP,DLP]=butter(LPorder,2*LP/fsample,'low');
DCO2f=filtfilt(NLP,DLP,DCO2);
DO2f=filtfilt(NLP,DLP,DO2);
end
if ~isempty(HP) && ~isempty(LP)
[NHP,DHP]=butter(HPorder,2*HP/fsample,'high');
[NLP,DLP]=butter(LPorder,2*LP/fsample,'low');
DCO2f=filtfilt(NLP,DLP,filtfilt(NHP,DHP,DCO2));
DO2f=filtfilt(NLP,DLP,filtfilt(NHP,DHP,DO2));
end
if isempty(HP) && isempty(LP)
DCO2f=DCO2;
DO2f=DO2;    
end

PETCO2filt=DCO2f;
PETO2filt=DO2f;
%% filter ASL and BOLD timecourses
DIM=size(ASL.CBF);
CBFfilt=NaN(DIM);
BOLDfilt=NaN(DIM);


tic
for ii=1:DIM(1)
  
    for jj=1:DIM(2)
        for kk=1:DIM(3)
            if ~isnan(ASL.CBF0(ii,jj,kk))  % run for all voxels within the mask
                
                cbf=squeeze(ASL.CBF(ii,jj,kk,:)); %CBF signal within the voxel
                bold=squeeze(BOLD.BOLDdet(ii,jj,kk,:)); %BOLD fractional and detrended change within the voxel

                    
                    %  filter fMRI data 
                    
                    if ~isempty(HP) && isempty(LP)              
                        Yf1=filtfilt(NHP,DHP,cbf);
                    end
                    if ~isempty(LP) && isempty(HP)       
                        Yf1=filtfilt(NLP,DLP,cbf);
                    end
                    if ~isempty(HP) && ~isempty(LP)
                        Yf1=filtfilt(NLP,DLP,filtfilt(NHP,DHP,cbf));
                    end
                    if isempty(HP) && isempty(LP)
                        Yf1=cbf;
                    end  
                    CBFfilt(ii,jj,kk,:)=Yf1-mean(Yf1);
                    bold(isnan(bold))=0;          
                    
                    Yb1=filtfilt(NHP,DHP,bold);
                    if ~isempty(HP) && isempty(LP)
                        Yb1=filtfilt(NHP,DHP,bold);
                    end
                    if ~isempty(LP) && isempty(HP)
                        Yb1=filtfilt(NLP,DLP,bold);
                    end
                    if ~isempty(HP) && ~isempty(LP)
                        Yb1=filtfilt(NLP,DLP,filtfilt(NHP,DHP,bold));
                    end
                    if isempty(HP) && isempty(LP)
                        Yb1=bold;
                    end
                    BOLDfilt(ii,jj,kk,:)=Yb1-mean(Yb1);                   
            end
        end
    end    
    disp(['Applying Temporal Filters to fMRI data...',num2str(round(100*ii/DIM(1))),'%']);
end
toc

%% Extract Median GM Signal
CBFGMfilt=[];
BOLDGMfilt=[];
  
    for j=1:size(CBFfilt,4)
        temp=squeeze(CBFfilt(:,:,:,j));
        temp1=temp(:);
        CBFGMfilt(j)=nanmedian(temp1(logical(ASL.GM(:))));
        temp=squeeze(BOLDfilt(:,:,:,j));
        temp1=temp(:);
        BOLDGMfilt(j)=nanmedian(temp1(logical(ASL.GM(:))));
    end
CBF0GM=nanmedian(ASL.CBF0(ASL.GM));

%% save data in the output structure

fMRI_filt.CBF0=ASL.CBF0;
fMRI_filt.GM=ASL.GM;
fMRI_filt.Parameters=Parameters;
fMRI_filt.PaO20=PaO20;
fMRI_filt.PaCO20=PaCO20;
fMRI_filt.PETO2filt=PETO2filt;
fMRI_filt.PETCO2filt=PETCO2filt;
fMRI_filt.CBFfilt=CBFfilt;
fMRI_filt.BOLDfilt=BOLDfilt;
fMRI_filt.CBFGMfilt=CBFGMfilt;
fMRI_filt.BOLDGMfilt=BOLDGMfilt;
fMRI_filt.CBF0GM=CBF0GM;

if FigOn==1
     FigTitle=Parameters.FigTitle;
    fs=10;

    FigH = figure('Position', get(0, 'Screensize'));
    subplot(2,2,1)
    plot([1:length(PETCO2)]/fsample,PETCO2)
    xlabel('Time (s)')
    ylabel('End-Tidal CO_{2} (mmHg)')
    set(gca,'FontSize',fs,'FontWeight','bold')
    set(gca,'FontSize',fs,'FontWeight','bold')
    title(FigTitle)
    subplot(2,2,2)
    plot([1:length(PETCO2)]/fsample,PETO2)
    xlabel('Time (s)')
    ylabel('End-Tidal O_{2} (mmHg)')
    set(gca,'FontSize',fs,'FontWeight','bold')
    subplot(2,2,3)
    plot([1:length(PETCO2)]/fsample,PETO2filt)
    hold on
    xlabel('Time (s)')
    ylabel('End-Tidal CO_{2} (mmHg)')
    set(gca,'FontSize',fs,'FontWeight','bold')
    title('Filtered')
    subplot(2,2,4)
    plot([1:length(PETCO2)]/fsample,PETO2filt)
    hold on
    xlabel('Time (s)')
    ylabel('End-Tidal O_{2} (mmHg)')
    set(gca,'FontSize',fs,'FontWeight','bold')
    title('Filtered')
    set(gcf,'color','w');
    set(gca,'FontSize',fs,'FontWeight','bold')

    FigH = figure('Position', get(0, 'Screensize'));
    subplot(1,2,1)
    plot([1:length(PETCO2)]/fsample,100*CBFGMfilt/CBF0GM)
    hold on
    xlabel('Time (s)')
    ylabel('CBF/CFB0 (%)')
    title([FigTitle,'; Filtered, GM Average'])
    set(gcf,'color','w');
    set(gca,'FontSize',fs,'FontWeight','bold')
    subplot(1,2,2)
    plot([1:length(PETCO2)]/fsample,100*BOLDGMfilt)
    hold on
    xlabel('Time (s)')
    ylabel('BOLD/BOLD0 (%)')
    set(gcf,'color','w');
    set(gca,'FontSize',fs,'FontWeight','bold')
    title('Filtered, GM Average')
end

    
  
        