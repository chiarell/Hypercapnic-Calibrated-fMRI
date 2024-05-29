
% Written by Chiarelli Atonio Maria, University G. D'Annunzio of Chieti-Pescara
% Original Date: 05/06/2021
%Last Update: 05/03/2024

%If used for publication, please site the following papers:

%Chiarelli, A. M., Germuska, M., Chandler, H., Stickland, R., Patitucci, E., Biondetti, E., ... & Wise, R. G. (2022). A flow-diffusion model of oxygen transport for quantitative mapping of cerebral metabolic rate of oxygen (CMRO2) with single gas calibrated fMRI. Journal of Cerebral Blood Flow & Metabolism, 42(7), 1192-1209.
%Germuska, M., Chandler, H. L., Stickland, R. C., Foster, C., Fasano, F., Okell, T. W., ... & Wise, R. G. (2019). Dual-calibrated fMRI measurement of absolute cerebral metabolic rate of oxygen consumption and effective oxygen diffusivity. Neuroimage, 184, 717-728.
%...

%The function inverts the Davis model extended with the flow-diffusion
%model of oxygen transport to estimate OEF0 and CMRO2 with a single
%calibrated fMRI experiment using a vasodilatatory signal (e.g. hypercapnia induced by CO2 inhalation or breath holding, or resting state data exploiting endogenous vascular signals). 
%The method used the General Linear Model (GLM) to estimate BOLD and ASL
%CVRs, the BOLD and ASL modulations and the BOLD and ASL SNRs
% The method  inverts the model and obtain OEF0,
% assuming a fixed a multiplicative constant equal to Aro/k of unit
% s^-1g^-beta dL^beta/(micromol/mmHg/ml/min) and a fixed mithocondrial
% pressure PmO2.

%Inputs:
%fMRIfilt structure output of the FilterCBF_BOLD_ET.m function
%GM Only: integer 1: compute values only within the GM mask, 0: compute
%values in the whole image (expect where NaNs are);
% RegressorType: vasodilalatory signal to be used as a regressor integer 1:
% 1: Filtered PetCO2; 2: z-scored filtered average GM BOLD signal (suggested for RS data); 3: z-scored filtered average GM CBF signal

%Parameters structure must report all the information in the example, with the
%same units
%     Parameters.TE=30;  %Time of echo of the BOLD sequence, units of ms
%     Parameters.TR=4.4; %Time of echo of the BOLD and ASL signal, units of s
%     Parameters.THP=150; %high pass cut-off time, units of s
%     Parameters.LHP=10; %low pass cut-off time, units of s
%     Parameters.beta=1.3; % adimensional
%     Parameters.alpha=0.38; % Grubb exponent
%     Parameters.phi=1.32; % Oxygen binding macity of hemoglobin, ml of oxygen/g of Hb
%     Parameters.epsilon=0.0031; % Oxygen plasma solubility, ml/mmHg/dL
%     Parameters.h=2.8; % Hill coefficient, adimensional
%     Parameters.Hb=14;  %hemoglobin in blood g/dL
%     OEFmin=0; % minimum OEF allowed
%     OEFmax=1; % maximum OEF allowed
%     OEFstep=10^-3; %OEF step for explicit search during the inversion (smaller steps slow the procedure, coarser steps make the estimate coarser)
%     Parameters.PmO2=0; % assumed average partial pressure at the mithocondria (mmHg)
%     %(units of mmHg, close to zero in Healthy subjects, i would use 10 in disease);
%     Parameters.Aro_k=8.8; % units of s^-1g^-beta dL^beta/(micromol/mmHg/ml/min) (suggested 10 s^-1g^-beta dL^beta/(micromol/mmHg/ml/min) at 3 T)
%     Parameters.n=10^5; %flow/cmro2 modulation ration. with isometabolism
%     n is very high.
%   Parameters.Mmin=0.01; % minimum M allowed
%   Parameters.Mmax=0.5; % maximum M allowed
%   Parameters.PmO2=0; % assumed average partial pressure at the mithocondria (mmHg)
%     %(units of mmHg, close to zero in Healthy subjects, i would use 10 in disease);
%   Parameters.Aro_k=8.8; % units of s^-1g^-beta dL^beta/(micromol/mmHg/ml/min) (suggested 10 s^-1g^-beta dL^beta/(micromol/mmHg/ml/min) at 3 T)
%   Parameters.n=10^5; %relative change in CBF over relative change in CMRO2
%     Parameters.shift=0;;% temporal (positive and negative) admissable time shift
%     (seconds)
%     between regressor and voxel signal allowed
% Parameter.InvMethod, string with either 'Search' (fminsearch matlab
% function used), 'Newton' (fminunc matlab function used), or 'Grid'
% (explicit grid search) (Initial guess for oef0=0.4)
%FigOn=1 creates matlab figures of the analysis outcome

%Outputs:
% The function returns voxelwise maps of OEF0 (Relative Units), CMR02 (micromol/100g/min), 
%CVR (%CBF/mmHg of CO2), BOLD_CVR (%BOLD/mmHg of CO2), 
%SNR_ASL of the hypermnic modulation and SNR_BOLD of the hypermnic modulation






function [outcome]=vasfMRI(fMRIfilt,GMonly,RegressorType,Parameters,FigOn)





% inizialize the outcome structure
DIM=size(fMRIfilt.CBFfilt);
outcome.M=NaN([DIM(1:3)]);
outcome.OEF0=NaN([DIM(1:3)]);
outcome.MCTT0=NaN([DIM(1:3)]);
outcome.CMRO20=NaN([DIM(1:3)]);
outcome.CVR=NaN([DIM(1:3)]);
outcome.BOLDCVR=NaN([DIM(1:3)]);
outcome.CVRlag=NaN([DIM(1:3)]);
outcome.BOLDCVRlag=NaN([DIM(1:3)]);
outcome.CVR_SNR=NaN([DIM(1:3)]);
outcome.BOLDCVR_SNR=NaN([DIM(1:3)]);
outcome.Parameters=Parameters;

% import the needed parameters
TR=Parameters.TR;
TE=Parameters.TE;

beta=Parameters.beta;
alpha=Parameters.alpha;
phi= Parameters.phi;
epsilon=Parameters.epsilon;
h=Parameters.h;
OEFmin=Parameters.OEFmin;
OEFmax=Parameters.OEFmax;
Mmin=Parameters.Mmin;
Mmax=Parameters.Mmax;
OEFstep=Parameters.OEFstep;
PmO2=Parameters.PmO2;
Aro_k=Parameters.Aro_k;
Hb=Parameters.Hb;
n=Parameters.n;
shift=Parameters.shift;
InvMethod=Parameters.InvMethod;

%% compute P50 from baseline CO2
PaCO20=fMRIfilt.PaCO20;
HCO3m=24; %mmol/L;
pH=6.1+log10(HCO3m/(0.03*PaCO20));
P50=221.87-26.37*pH; % blood oxygen tension at which haemoglobin is 50% saturated (mmHg)

%% compute saturation and oxygen concentration in arterial blood at rest
PaO20=fMRIfilt.PaO20;
SaO20=1/(1+(P50/PaO20)^h);
CaO20=phi*SaO20*Hb+epsilon*PaO20; % ml of oxygen/dl of blood
PaO2m=PaO20; % mmHg %%%%%%%%%% m stands for modulation value, WARNING: Modify if you assume changes in PaO20 with hypercapnia
SaO2m=1/(1+(P50/PaO2m)^h);
CaO2m=phi*SaO2m*Hb+epsilon*PaO2m; % ml of oxygen/dl of blood

outcome.PaCO20=PaCO20;
outcome.PaO20=PaO20;
outcome.CaO20=CaO20;
outcome.SaO20=SaO20;
outcome.PaO2m=PaO2m; 
outcome.SaO2m=SaO2m;
outcome.CaO2m=CaO2m; 
outcome.P50=P50;

%% select regressor

CBF0=fMRIfilt.CBF0;
CBFfilt=fMRIfilt.CBFfilt;
BOLDfilt=fMRIfilt.BOLDfilt;
if GMonly==1
    MASK=fMRIfilt.GM;
else
    MASK=~isnan(CBF0);
end

if RegressorType==1
    regressor=fMRIfilt.PETCO2filt;
end
if RegressorType==2
    regressor=zscore(fMRIfilt.BOLDGMfilt);
end
if RegressorType==3
    regressor=zscore(fMRIfilt.CBFGMfilt);
end
fsample=1/TR;
shiftS=round(shift*fsample);

%% compute CVR, BOLD-CVR and fit the flow diffusion model to compute OEF0 and CMRO20 (Chiarelli et al., JCBFM, 2022)

OEF0Grid=OEFmin:OEFstep:OEFmax;
DIM=size(fMRIfilt.CBFfilt);
tic 
for ii=1:DIM(1)

    for jj=1:DIM(2)
        for kk=1:DIM(3)
            if MASK(ii,jj,kk)>0  % run for all voxels within the mask


                if shiftS==0
                    [CBFb,~,statsF]=glmfit(regressor,squeeze(CBFfilt(ii,jj,kk,:)));
                    [BOLDb,~,statsB]=glmfit(regressor,squeeze(BOLDfilt(ii,jj,kk,:)));
                    indF=0;
                    indB=0;
                else
                    [z,lags]=xcorr(regressor,squeeze(CBFfilt(ii,jj,kk,:)),shiftS);
                    indF=lags(z==max(z));
                    [CBFb,~,statsF]=glmfit(regressor,circshift(squeeze(CBFfilt(ii,jj,kk,:)),indF(1)));
                    [z,lags]=xcorr(regressor,squeeze(BOLDfilt(ii,jj,kk,:)),shiftS);
                    indB=lags(z==max(z));
                    [BOLDb,~,statsB]=glmfit(regressor,circshift(squeeze(BOLDfilt(ii,jj,kk,:)),indB(1)));
                end



                % compute SNR based on beta-weights of the GLM
                outcome.CVR_SNR(ii,jj,kk)=statsF.t(2);
                outcome.BOLDCVR_SNR(ii,jj,kk)=statsB.t(2);
                % invert the flow diff model
                BOLDm=(BOLDb(2))*(max(regressor));
                CBF0i=CBF0(ii,jj,kk);
                
                CBFm=CBF0i+CBFb(2)*(max(regressor));
                
                [Mest, OEF0est, MCTT0est,CMRO20est]=Inverse_BOLD_OxFlowDiff(OEF0Grid,Mmin,Mmax,BOLDm,CBFm,CBF0i,CaO2m,CaO20,Hb,n,phi,h,alpha,beta,TE,P50,Aro_k,PmO2,InvMethod);

                outcome.M(ii,jj,kk)=Mest;
                outcome.MCTT0(ii,jj,kk)=MCTT0est;
                outcome.MCTT0(ii,jj,kk)=MCTT0est;
                outcome.OEF0(ii,jj,kk)=OEF0est;
                outcome.CMRO20(ii,jj,kk)=CMRO20est;

                outcome.CVR(ii,jj,kk)=(CBFm/CBF0i-1)/max(regressor);
                outcome.BOLDCVR(ii,jj,kk)=BOLDb(2);
                outcome.CBF0(ii,jj,kk)=CBF0i;
                outcome.CVRlag(ii,jj,kk)=-indF(1)/fsample;
                outcome.BOLDCVRlag(ii,jj,kk)=-indB(1)/fsample;

             
            end
        end
    end 
disp(['Computing Maps...',num2str(round(100*ii/DIM(1))),'%']);
end
toc
%% save other information in the new structure
outcome.CBF0=CBF0;
outcome.regressor=regressor;
outcome.RegressorType=RegressorType;
outcome.Parameters=Parameters;
outcome.GM=fMRIfilt.GM;
outcome.CBF0GM=nanmedian(outcome.CBF0(outcome.GM));
outcome.MGM=nanmedian(outcome.M(outcome.GM));
outcome.OEF0GM=nanmedian(outcome.OEF0(outcome.GM));
outcome.CMRO20GM=nanmedian(outcome.CMRO20(outcome.GM));
outcome.MCTT0GM=nanmedian(outcome.MCTT0(outcome.GM));
outcome.BOLDCVRGM=nanmedian(outcome.BOLDCVR(outcome.GM));
outcome.CVRGM=nanmedian(outcome.CVR(outcome.GM));
outcome.BOLDCVR_SNRGM=nanmedian(outcome.BOLDCVR_SNR(outcome.GM));
outcome.CVR_SNRGM=nanmedian(outcome.CVR_SNR(outcome.GM));

if FigOn==1
    FigTitle=Parameters.FigTitle;
    
    fs=20;
    DIM=size(outcome.M);
     
    FigH = figure('Position', get(0, 'Screensize'));
    subplot(2,2,1)
    imagesc(100*rot90(outcome.BOLDCVR(:,:,round(DIM(3)/2))))
    axis off
    axis equal
    caxis([0 100*prctile(outcome.BOLDCVR(:),95)])
    title([FigTitle,'; BOLD CVR (A.U.)'])
    hcb=colorbar;
    title(hcb,'(%)')
    set(gca,'FontSize',fs,'FontWeight','bold')
    subplot(2,2,2)
    A=100*outcome.BOLDCVR.*~isnan(outcome.OEF0);
    A(A>prctile(A,95))=NaN;
    A(A<=0)=NaN;
    hist(A(:),round(sqrt(sum(~isnan(A(:))))))
    axis([0 100*prctile(outcome.BOLDCVR(:),95) 0 Inf])
    xlabel('BOLD CVR (A.U.)')
    ylabel('Frequency')
    set(gca,'FontSize',fs,'FontWeight','bold') 
    axis square
    subplot(2,2,3)
    imagesc(100*rot90(outcome.CVR(:,:,round(DIM(3)/2))))
    axis off
    axis equal
    caxis([0 100*prctile(outcome.CVR(:),95)])
    title('CVR (A.U.)')
    hcb=colorbar;
    title(hcb,'(%)')
    set(gca,'FontSize',fs,'FontWeight','bold')
    subplot(2,2,4)
    A=100*outcome.CVR.*~isnan(outcome.OEF0);
    A(A>prctile(A,95))=NaN;
    A(A<=0)=NaN;
    hist(A(:),round(sqrt(sum(~isnan(A(:))))))
    axis([0 100*prctile(outcome.CVR(:),95) 0 Inf])
    xlabel('CVR (A.U.)')
    ylabel('Frequency')
    set(gca,'FontSize',fs,'FontWeight','bold') 
    set(gcf,'color','w'); 
    axis square



    FigH = figure('Position', get(0, 'Screensize'));
    subplot(3,2,1)
    imagesc(100*rot90(outcome.M(:,:,round(DIM(3)/2))))
    axis off
    axis equal
    caxis([0 100*prctile(outcome.M(:),95)])
    title('M')
    hcb=colorbar;
    title(hcb,'(%)')
    set(gca,'FontSize',fs,'FontWeight','bold')
    subplot(3,2,2)
    A=100*outcome.M;
    hist(A(:),round(sqrt(sum(~isnan(A(:))))))
    axis([0 100*prctile(outcome.M(:),95) 0 Inf])
    xlabel('M (%)')
    ylabel('Frequency')
    set(gca,'FontSize',fs,'FontWeight','bold') 
    axis square
    title(FigTitle)
    subplot(3,2,3)
    imagesc(rot90(outcome.OEF0(:,:,round(DIM(3)/2))))
    axis off
    axis equal
    caxis([0 1])
    hcb=colorbar;
    title('OEF_{0}')
    set(gca,'FontSize',fs,'FontWeight','bold')
    subplot(3,2,4)
    A=outcome.OEF0;
    hist(A(:),round(sqrt(sum(~isnan(A(:))))))
    axis([0 1 0 Inf])
    xlabel('OEF_{0}')
    ylabel('Frequency')
    set(gca,'FontSize',fs,'FontWeight','bold') 
    axis square
    subplot(3,2,5)
    imagesc(rot90(outcome.CMRO20(:,:,round(DIM(3)/2))))
    axis off
    axis equal
    caxis([0 prctile(outcome.CMRO20(:),95)])
    hcb=colorbar;
    title('CMRO_{2,0}')
    title(hcb,'(\mumol/100g/min)')
    set(gca,'FontSize',fs,'FontWeight','bold')
    subplot(3,2,6)
    A=outcome.CMRO20;
    %A(A>300)=NaN;
    hist(A(:),round(sqrt(sum(~isnan(A(:))))))
    axis([0 prctile(outcome.CMRO20(:),95) 0 Inf])
    xlabel('CMRO_{2,0} (\mumol/100g/min)')
    ylabel('Frequency')
    set(gca,'FontSize',fs,'FontWeight','bold') 
    set(gcf,'color','w'); 
    axis square
     
     
     
     
     
     
     
     
     




end

end
  
%% Auxillary functions
% [BOLD, M, MCTT]=BOLD_OxFlowDiff(CBFm,CBF0,CaO2m,CaO20,OEF0,Hb,n,phi,h,alpha,beta,TE,P50,Aro_k,PmO2)
% Written by Chiarelli Atonio Maria, University G. D'Annunzio of Chieti-Pescara
% Original Date: 18/03/2024
% Last UpDate: 18/03/2024

%The function computes the forward model of the BOLD signal integrated with
%the flow-diffusion model of oxygen transport
%Inputs:
%CBFm CBF during modulation in ml/100g/min
%CBF0 CBF at baseline in ml/100g/min
%CaO2m oxygen concentration in blood during modulation in ml of oxygen/dl of blood
% CaO20 oxygen concentration at baseline in ml of oxygen/dl of blood
% OEF0 baseline oxygen extraction fraction
%hemoglobin in blood g/dL
%n relative change in CBF over relative change in CMRO2 (put a very high
%value if isometabolism is assumed, e.g. 10^5).
% phi=1.32; % Oxygen binding macity of hemoglobin, ml of oxygen/g of Hb
%h=2.8; % Hill coefficient, adimensional
%alpha=0.38 Grubb constant
%beta=1.3 % beta od the davis model
% TE echo time in ms
% P50, oxygen pressure in blood at 50% hemoglobin saturation 
%Aro_k=8.8; % units of s^-1g^-beta dL^beta/(micromol/mmHg/ml/min)
% PmO2 oxygen pressure at the mitochondria
% Output: BOLD signal, maximum BOLD modulation, Mean capillary transit time
% in seconds  
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
