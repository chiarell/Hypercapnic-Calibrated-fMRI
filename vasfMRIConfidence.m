
% Written by Chiarelli Atonio Maria, University G. D'Annunzio of Chieti-Pescara
% Original Date: 05/01/2024
%Last Update: 05/03/2024



%The function uses experimental noise within a simulation framework to
%understand the confidence in parameters estimation from the experimental
%BOLD-ASL acquisitions
%Inputs:
%outcome: structure output of hcfmri.m
%fMRIfilt structure output of the FilterCBF_BOLD_ET.m function
%ASL: structure output of RAW2BF.m
%n_iterations: number of iterations within the simulation framework
%Fig on: 1 create figures, 0 do not create figures



function [outcome]=vasfMRIConfidence(outcome,fMRI_filt,ASL,n_iterations,FigOn)


%% import the needed parameters
THP=fMRI_filt.Parameters.THP;
TLP=fMRI_filt.Parameters.TLP;
HPorder=fMRI_filt.Parameters.HPorder;
LPorder=fMRI_filt.Parameters.LPorder;

HP=1./THP;
LP=1./TLP;

% import the needed parameters
TR=outcome.Parameters.TR;
TE=outcome.Parameters.TE;

beta=outcome.Parameters.beta;
alpha=outcome.Parameters.alpha;
phi=outcome.Parameters.phi;
h=outcome.Parameters.h;
OEFmin=outcome.Parameters.OEFmin;
OEFmax=outcome.Parameters.OEFmax;
Mmin=outcome.Parameters.Mmin;
Mmax=outcome.Parameters.Mmax;
OEFstep=outcome.Parameters.OEFstep;
PmO2=outcome.Parameters.PmO2;
Aro_k=outcome.Parameters.Aro_k;
Hb=outcome.Parameters.Hb;
n=outcome.Parameters.n;
shift=outcome.Parameters.shift;
InvMethod=outcome.Parameters.InvMethod;
PaCO20=outcome.PaCO20;
CaO20=outcome.CaO20;
CaO2m=CaO20;


fsample=1/TR;
shiftS=round(shift*fsample);
DIMt=size(fMRI_filt.CBFfilt,4);
CBF0GM=fMRI_filt.CBF0(fMRI_filt.GM);
CBF0_SNRGM=(ASL.tSNR(fMRI_filt.GM))*sqrt(length(fMRI_filt.Parameters.baseline));
CVRGM=outcome.CVR(fMRI_filt.GM);
CVR_SNRGM=outcome.CVR_SNR(fMRI_filt.GM);
BOLDCVR_SNRGM=outcome.BOLDCVR_SNR(fMRI_filt.GM);




%% define filter
if ~isempty(HP) && isempty(LP)
[NHP,DHP]=butter(HPorder,2*HP/fsample,'high');
end
if ~isempty(LP) && isempty(HP)
[NLP,DLP]=butter(LPorder,2*LP/fsample,'low');

end
if ~isempty(HP) && ~isempty(LP)
[NHP,DHP]=butter(HPorder,2*HP/fsample,'high');
[NLP,DLP]=butter(LPorder,2*LP/fsample,'low');
end
%%

%% compute P50 from baseline CO2
HCO3m=24; %mmol/L;
pH=6.1+log10(HCO3m/(0.03*PaCO20));
P50=221.87-26.37*pH; % blood oxygen tension at which haemoglobin is 50% saturated (mmHg)


%% understand true SNRs

% for q=1:50
%     q
% CVR_SNRtry(:,q)=CVR_SNRGM-5+q*0.1;
% BOLDCVR_SNRtry(:,q)=BOLDCVR_SNRGM-5+q*0.1;
% for i=1:10^3
% regressorRaw=randn(DIMt,1);
% CBFRaw=randn(DIMt,1);
% BOLDRaw=randn(DIMt,1);
% 
% ind=randperm(size(CVR_SNRtry,1));
% CVRtry=CVRGM(ind(1));
% if ~isempty(HP) && isempty(LP)
%     regressor=filtfilt(NLP,DLP,filtfilt(NHP,DHP,regressorRaw));
%     CBFm=filtfilt(NLP,DLP,filtfilt(NHP,DHP,CBFRaw));
%     BOLDm=filtfilt(NLP,DLP,filtfilt(NHP,DHP,BOLDRaw));
% end
% if ~isempty(LP) && isempty(HP)
%     regressor=filtfilt(NLP,DLP,regressorRaw);
%     CBFm=filtfilt(NLP,DLP,CBFRaw);
%     BOLDm=filtfilt(NLP,DLP,BOLDRaw);
% end
% if ~isempty(HP) && ~isempty(LP)
%     regressor=filtfilt(NHP,DHP,regressorRaw);
%     CBFm=filtfilt(NLP,DLP,CBFRaw);
%     BOLDm=filtfilt(NLP,DLP,BOLDRaw);
% end
% if isempty(HP) && isempty(LP)
%     regressor=regressorRaw;
%     CBFm=CBFRaw;
%     CBFm=BOLDRaw;
% end
% regressor=regressor/std(regressor);
% CVRr=CVRtry*regressor;
% CBFsignal=CVRr+CBFm/CVR_SNRtry(ind(1),q)*sqrt(size(ASL.CBF,4))*std(CVRr)/std(CBFm);
% BOLDsignal=CVRr+BOLDm/BOLDCVR_SNRtry(i)*sqrt(size(ASL.CBF,4))*std(CVRr)/std(BOLDm);
% 
%     [z,lags]=xcorr(regressor,CBFsignal,shiftS);
%     indF=lags(z==max(z));
%     [b,m,stats]=glmfit(regressor,circshift(CBFsignal,indF(1)));
%     CVR_SNRe(i)=stats.t(2);
% 
%     [z,lags]=xcorr(regressor,BOLDsignal,shiftS);
%     indF=lags(z==max(z));
%     [b,m,stats]=glmfit(regressor,circshift(BOLDsignal,indF(1)));
%     BOLDCVR_SNRe(i)=stats.t(2);
% end
% biasCVR(q)=nanmean(CVR_SNRe(:,q))-nanmean(CVR_SNRGM);
% biasBOLDCVR(q)=nanmean(BOLDCVR_SNRe(:,q))-nanmean(BOLDCVR_SNRGM);
% end
% 
% [m ind]=find(abs(biasCVR)==min(abs(biasCVR)))
% CVR_SNRGM=CVR_SNRGM-5+ind*0.1;
% [m ind]=find(abs(biasBOLDCVR)==min(abs(biasBOLDCVR)))
% BOLDCVR_SNRGM=BOLDCVR_SNRGM-5+ind*0.1;

%%

for i=1:n_iterations

M(i)=Mmin-1;
while M(i)<Mmin || M(i)>Mmax

OEF0(i)=0;
while(OEF0(i)<0.1) || (OEF0(i)>0.9)
    OEF0(i)=0.4+randn(1)/8;
end

CVR(i)=-1;
while CVR(i)<0 || CVR(i)>0.5
ind=randperm(length(CVRGM));
CVR(i)=CVRGM(ind(1));
end
CBF0(i)=CBF0GM(ind(1));
CMRO20(i)=CBF0(i)*OEF0(i)*CaO20*10^-2*10^-3*10^6/22.4; % compute CMRO2 in micromol/100g/min
CBF0_SNR(i)=CBF0_SNRGM(ind(1));

CVR_SNR(i)=CVR_SNRGM(ind(1));
BOLDCVR_SNR(i)=BOLDCVR_SNRGM(ind(1));

CBF0est(i)=CBF0(i)+CBF0(i)*randn(1)/CBF0_SNR(i);


regressorRaw=randn(DIMt,1);
CBFRaw=randn(DIMt,1);
CBFNCRaw=randn(DIMt,1);
BOLDRaw=randn(DIMt,1);

if ~isempty(HP) && isempty(LP)
    regressor=filtfilt(NLP,DLP,filtfilt(NHP,DHP,regressorRaw));
    CBFm=filtfilt(NLP,DLP,filtfilt(NHP,DHP,CBFRaw));
 
    BOLDm=filtfilt(NLP,DLP,filtfilt(NHP,DHP,BOLDRaw));
end
if ~isempty(LP) && isempty(HP)
    regressor=filtfilt(NLP,DLP,regressorRaw);
    CBFm=filtfilt(NLP,DLP,CBFRaw);
  
    BOLDm=filtfilt(NLP,DLP,BOLDRaw);
end
if ~isempty(HP) && ~isempty(LP)
    regressor=filtfilt(NHP,DHP,regressorRaw);
    CBFm=filtfilt(NLP,DLP,CBFRaw);

    BOLDm=filtfilt(NLP,DLP,BOLDRaw);
end
if isempty(HP) && isempty(LP)
    regressor=regressorRaw;
    CBFm=CBFRaw;
    BOLDm=BOLDRaw;
end

% CBFm=CBFfiltGM(ind(1),:)';
% BOLDm=BOLDfiltGM(ind(1),:)';

regressor=regressor/std(regressor);



CVRr=CVR(i)*regressor;
CBFsignal=CVRr+CBFm/CVR_SNR(i)*sqrt(size(ASL.CBF,4))*std(CVRr)/std(CBFm);


if shiftS==0
    [b,m,stats]=glmfit(regressor,CBFsignal);
else
    [z,lags]=xcorr(regressor,CBFsignal,shiftS);
    indF=lags(z==max(z));
    [b,m,stats]=glmfit(regressor,circshift(CBFsignal,indF(1)));
   
end
CVRest(i)=b(2);
CVRlag(i)=indF(1);
CVR_SNRest(i)=stats.t(2);



[BOLDCVR(i),M(i), MCTT0(i)]=BOLD_OxFlowDiff(CBF0(i)*(1+CVR(i)),CBF0(i),CaO2m,CaO20,OEF0(i),Hb,n,phi,h,alpha,beta,TE,P50,Aro_k,PmO2);

end

BOLDsignal=squeeze(BOLDCVR(i))*regressor+BOLDm/BOLDCVR_SNR(i)*sqrt(size(ASL.CBF,4))*std(BOLDCVR(i)*regressor)/std(BOLDm);

if shiftS==0
    [b,m,stats]=glmfit(regressor,BOLDsignal);
else
    [z,lags]=xcorr(regressor,BOLDsignal,shiftS);
    indB=lags(z==max(z));
    [b,m,stats]=glmfit(regressor,circshift(BOLDsignal,indB(1)));
end
BOLDCVRest(i)=b(2);
BOLDCVRlag(i)=indB(1);
BOLDCVR_SNRest(i)=stats.t(2);

OEF0Grid=OEFmin:OEFstep:OEFmax;

if BOLDCVR(i)<0.1 && BOLDCVRest(i)>0 && BOLDCVRest(i)<0.1 && CVRest(i)<0.5
[Mest(i), OEF0est(i), MCTT0est(i),CMRO20est(i)]=Inverse_BOLD_OxFlowDiff(OEF0Grid,Mmin,Mmax,BOLDCVRest(i),CBF0est(i)*(1+CVRest(i)),CBF0est(i),CaO2m,CaO20,Hb,n,phi,h,alpha,beta,TE,P50,Aro_k,PmO2,InvMethod);
  Mest1(i)=BOLDCVRest(i)./CVRest(i); 

else
    CVR(i)=NaN;
BOLDCVR(i)=NaN;
BOLDCVRest(i)=NaN;
CVRest(i)=NaN;
Mest(i)=NaN;
OEF0est(i)=NaN;
MCTT0est(i)=NaN;
CMRO20est(i)=NaN;
end
if rem(100*i/n_iterations,1)==0
disp(['Simulating data with experimental SNR...',num2str(round(100*i/n_iterations)),'%']);
end

regressors(:,i)=regressor;
CBFsignals(:,i)=CBFsignal;
BOLDsignals(:,i)=BOLDsignal;
end

outcome.Simulations.regressors=regressors;
outcome.Simulations.CBFsignals=CBFsignals;
outcome.Simulations.BOLDsignals=BOLDsignals;
outcome.Simulations.BOLDCVR=BOLDCVR;
outcome.Simulations.BOLDCVRest=BOLDCVRest;
outcome.Simulations.BOLDCVR_SNR=BOLDCVR_SNR;
outcome.Simulations.BOLDCVR_SNRest=BOLDCVR_SNRest;
outcome.Simulations.CVR=CVR;
outcome.Simulations.CVRest=CVRest;
outcome.Simulations.CVR_SNR=CVR_SNR;
outcome.Simulations.CVR_SNRest=CVR_SNRest;
outcome.Simulations.CBF0=CBF0;
outcome.Simulations.CBF0est=CBF0est;
outcome.Simulations.CBF0_SNR=CBF0_SNR;
outcome.Simulations.M=M;
outcome.Simulations.Mest=Mest;
outcome.Simulations.OEF0=OEF0;
outcome.Simulations.OEF0est=OEF0est;
outcome.Simulations.CMRO20=CMRO20;
outcome.Simulations.CMRO20est=CMRO20est;

if FigOn==1
FigTitle=outcome.Parameters.FigTitle;
   FigH = figure('Position', get(0, 'Screensize'));
   subplot(2,3,1)
scatter(100*BOLDCVR,100*BOLDCVRest)%,'k','MarkerFaceColor', 'k')
hold on
plot([0 100], [0 100],'k','LineWidth',3)
RMSE=100*sqrt(nanmean((BOLDCVR-BOLDCVRest).^2))./nanmean(BOLDCVR);
title(['Voxel RMSE: ',num2str(RMSE),'%'])
xlabel('BOLD CVR (%, Simulated)')
ylabel('BOLD CVR (%, Estimated)')
set(gca,'FontSize',20,'FontWeight','bold')
set(gcf,'color','w');
h1=lsline;
h1.LineWidth=3;
h1.Color='r';
axis([0 1.05*max(100*BOLDCVR) 0 1.05*max(100*BOLDCVR)])
axis square
box on
  subplot(2,3,2)
scatter(100*CVR,100*CVRest)%,'k','MarkerFaceColor', 'k')
hold on
plot([0 100], [0 100],'k','LineWidth',3)
RMSE=100*sqrt(nanmean((CVR-CVRest).^2))./nanmean(CVR);
title(['Voxel RMSE: ',num2str(RMSE),'%'])
xlabel('CVR (%, Simulated)')
ylabel('CVR (%, Estimated)')
set(gca,'FontSize',20,'FontWeight','bold')
set(gcf,'color','w');
h1=lsline;
h1.LineWidth=3;
h1.Color='r';
axis([0 1.05*max(100*CVR) 0 1.05*max(100*CVR)])
axis square
box on
  subplot(2,3,3)
scatter(100*M,100*Mest)%,'k','MarkerFaceColor', 'k')
hold on
plot([0 100], [0 100],'k','LineWidth',3)
RMSE=100*sqrt(nanmean((M-Mest).^2))./nanmean(M);
title(['Voxel RMSE: ',num2str(RMSE),'%'])
xlabel('M (%, Simulated)')
ylabel('M (%, Estimated)')
set(gca,'FontSize',20,'FontWeight','bold')
set(gcf,'color','w');
h1=lsline;
h1.LineWidth=3;
h1.Color='r';
axis([0 1.05*max(100*M) 0 1.05*max(100*M)])
axis square
box on
 subplot(2,3,4)
scatter(100*OEF0,100*OEF0est)%,'k','MarkerFaceColor', 'k')
hold on
plot([0 100], [0 100],'k','LineWidth',3)
RMSE=100*sqrt(nanmean((OEF0-OEF0est).^2))./nanmean(OEF0);
title(['RMSE: ',num2str(RMSE),'%'])
xlabel('OEF_{0} (%, Simulated)')
ylabel('OEF_{0} (%, Estimated)')
set(gca,'FontSize',20,'FontWeight','bold')
set(gcf,'color','w');
h1=lsline;
h1.LineWidth=3;
h1.Color='r';
axis([0 1.05*max(100*OEF0) 0 1.05*max(100*OEF0)])
axis square
box on
subplot(2,3,5)
scatter(CBF0,CBF0est)%,'k','MarkerFaceColor', 'k')
hold on
plot([0 300], [0 300],'k','LineWidth',3)
RMSE=100*sqrt(nanmean((CBF0-CBF0est).^2))./nanmean(CBF0);
title(['Voxel RMSE: ',num2str(RMSE),'%'])
xlabel('CBF_{0} (ml/100g/min, Simulated)')
ylabel('CBF_{0} (mil/100g/min, Estimated)')
set(gca,'FontSize',20,'FontWeight','bold')
set(gcf,'color','w');
h1=lsline;
h1.LineWidth=3;
h1.Color='r';
axis([0 1.05*max(CBF0) 0 1.05*max(CBF0)])
axis square
box on
 subplot(2,3,6)
scatter(CMRO20,CMRO20est)%,'k','MarkerFaceColor', 'k')
hold on
plot([0 100], [0 100],'k','LineWidth',3)
RMSE=100*sqrt(nanmean((CMRO20-CMRO20est).^2))./nanmean(CMRO20);
title(['Voxel RMSE: ',num2str(RMSE),'%'])
xlabel('CMRO_{2,0} (\mumol/100g/min, Simulated)')
ylabel('CMRO_{2,0} (\mumol/100g/min, Estimated)')
set(gca,'FontSize',20,'FontWeight','bold')
set(gcf,'color','w');
h1=lsline;
h1.LineWidth=3;
h1.Color='r';
axis([0 1.05*max(CMRO20) 0 1.05*max(CMRO20)])
axis square
box on


end
end

%% Auxillary functions

function [BOLD, M, MCTT]=BOLD_OxFlowDiff(CBFm,CBF0,CaO2m,CaO20,OEF0,Hb,n,phi,h,alpha,beta,TE,P50,Aro_k,PmO2)


SvO20=CaO20./(phi*Hb).*(1-OEF0);

m=1+(CBFm/CBF0-1)/n;
OEFm=m*CBF0.*OEF0.*CaO20./(CBFm.*CaO2m);
SvO2m=CaO2m./(phi*Hb).*(1-OEFm);
Dm=(1-(CBFm./CBF0).^alpha.*((1-SvO2m)./(1-SvO20)).^beta);
Tm=OEF0*CaO20/((P50*(2/OEF0-1).^1/h)-PmO2);
M=1/(2.24*10^5)*Aro_k*TE*CBF0*Tm*((1-SvO20)*Hb).^beta; % compute Maximum bold modulation in relative units
MCTT=6*10^6/3*M/(Aro_k*TE*CBF0*((1-SvO20)*Hb).^beta);  % compute mean capillary transit time in seconds
BOLD=M*Dm;

end  
%  mse function to be used for the inverse of the flow diffusion modelling function

function [MSE]=MSE_BOLD_OxFlowDiff(OEF0,BOLDm,CBFm,CBF0,CaO2m,CaO20,Hb,n,phi,h,alpha,beta,TE,P50,Aro_k,PmO2)
[BOLDe]=BOLD_OxFlowDiff(CBFm,CBF0,CaO2m,CaO20,OEF0,Hb,n,phi,h,alpha,beta,TE,P50,Aro_k,PmO2);
MSE=(BOLDm-BOLDe)^2;
end
% inverse function
function [Mest, OEF0est, MCTT0est,CMRO20est]=Inverse_BOLD_OxFlowDiff(OEF0Grid,Mmin,Mmax,BOLDm,CBFm,CBF0,CaO2m,CaO20,Hb,n,phi,h,alpha,beta,TE,P50,Aro_k,PmO2,Method)


if strcmp(Method,'Search')
    fun = @(OEF0)MSE_BOLD_OxFlowDiff(OEF0,BOLDm,CBFm,CBF0,CaO2m,CaO20,Hb,n,phi,h,alpha,beta,TE,P50,Aro_k,PmO2);
    OEF0guess =0.4;
    OEF0est = fminsearch(fun,OEF0guess);
    [~, Mest, MCTT0est]=BOLD_OxFlowDiff(CBFm,CBF0,CaO2m,CaO20,OEF0est,Hb,n,phi,h,alpha,beta,TE,P50,Aro_k,PmO2);
    CMRO20est=OEF0est.*CBF0*CaO20*10^-2*10^-3*10^6/22.4; % compute CMRO2 in micromol/100g/min
end

if strcmp(Method,'Newton')
    fun = @(OEF0)MSE_BOLD_OxFlowDiff(OEF0,BOLDm,CBFm,CBF0,CaO2m,CaO20,Hb,n,phi,h,alpha,beta,TE,P50,Aro_k,PmO2);
    OEF0guess =0.4;
    OEF0est = fminunc(fun,OEF0guess);
    [~, Mest, MCTT0est]=BOLD_OxFlowDiff(CBFm,CBF0,CaO2m,CaO20,OEF0est,Hb,n,phi,h,alpha,beta,TE,P50,Aro_k,PmO2);
    CMRO20est=OEF0est.*CBF0*CaO20*10^-2*10^-3*10^6/22.4; % compute CMRO2 in micromol/100g/min
end

if strcmp(Method,'Grid')

    MSE=NaN(size(OEF0Grid));
    Me=NaN(size(OEF0Grid));
    MCTT0e=NaN(size(OEF0Grid));
    for i=1:length(OEF0Grid)
        [BOLDe, Me(i), MCTT0e(i)]=BOLD_OxFlowDiff(CBFm,CBF0,CaO2m,CaO20,OEF0Grid(i),Hb,n,phi,h,alpha,beta,TE,P50,Aro_k,PmO2);
        MSE(i)=(BOLDm-BOLDe)^2;
    end
    [~, OEF_ind]=min(MSE);
    Mest=real(Me(OEF_ind));
    MCTT0est=real(MCTT0e(OEF_ind));
    OEF0est=OEF0Grid(OEF_ind);
    CMRO20est=OEF0est.*CBF0*CaO20*10^-2*10^-3*10^6/22.4; % compute CMRO2 in micromol/100g/min
end

if ~strcmp(Method,'Search') && ~strcmp(Method,'Newton') && ~strcmp(Method,'Grid')
    disp('Inversion method unavailable, please choose between Search, Newton or Grid')
    Mest=NaN;
    MCTT0est=NaN;
    OEF0est=NaN;
    CMRO20est=NaN;
else
    if (OEF0est<=OEF0Grid(1)) || (OEF0est>=OEF0Grid(end)) || Mest<Mmin || Mest>Mmax
        Mest=NaN;
        MCTT0est=NaN;
        OEF0est=NaN;
        CMRO20est=NaN;
    end

end
end









