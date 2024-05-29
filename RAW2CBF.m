% Written by Chiarelli Antonio Maria, University G. D'Annunzio of Chieti-Pescara
% Original Date: 08/03/2024
% Last Update: 05/03/2024



%The function estimate CBF timecourses and average CBF map from from BOLD-ASL fMRI.
% The volumes should be preprocessed including susceptibility distortions
% correction, motion correction and registration to M0
%Inputs:
% PathIn ASL temporal volumes
% M0BrainPathIn: M0 brain masked map
% tag: 1 if first volume is tag, 0 if first volume is control
%Parameters.baseline, baseline points to evaluate CBF0
% Parameters structure with the following field
% Parameters.lambda, %brain/blood partition coefficient;
% Parameters.PLD, post label delay (s)
% Parameters.T1b,  T1 of blood (s)
%Parameters.alphaASL, tagging inversion efficiency
%Parameters.alphaASLInv, efficiency correction due to background suppression
%Parameters.tau, tagging time (s)
%Parameters.TR, Effective Repetition Time of dexi (seconds)
%Parameters.GMTH, an algorithm to extract the GM PVE directly from CBF is
%implemented, GM threhsold is used to define the GM Mask
%Parameters.FigTitle= String with the title of the figures
% FigOn, 1 plot figures, 0 do not plot figures
%Output:
% ASL structure with the following fields:

%ASL.CBF, 4D CBF (x,y,z,t)
%ASL.DM0, average perfusion signal withoout M0 normalization 
%ASL.PERFUSION0, average perfusion signal with M0 normalization 
%ASL.CBF0, average CBF
%ASL.CBFGM, GM Median CBF timecourse
%ASL.PVEGM, GM PVE
%ASL.GM, GM mask
%ASL.tSNR,; CBF temporal SNR
%ASL.Parameters, Parameters used for computation







function [ASL]=RAW2CBF(PathIn,M0BrainPathIn,tag,Parameters,FigOn)

%% Import Images
TE1_hdr = spm_vol(PathIn); 
TE1=spm_read_vols(TE1_hdr);

M0b_hdr = spm_vol(M0BrainPathIn);
M0b=spm_read_vols(M0b_hdr);

%% Import Parameters
baseline=Parameters.baseline;
lambda=Parameters.lambda; %brain/blood partition coefficient;
PLD= Parameters.PLD;% post label delay (s)
T1b=Parameters.T1b; % T1 of blood (s)
alphaASL= Parameters.alphaASL; % tagging inversion efficiency
alphaASLInv=Parameters.alphaASLInv; % efficiency correction due to background suppression
tau=Parameters.tau; % tagging time (s)
TR=Parameters.TR; % Effective Repetition Time of dexi (s)
GMTH=Parameters.GMTH;% GM threhsold


%% Perform surround subtraction on raw TE1
D=diff(TE1,1,4);
if tag==1
    D(:,:,:,2:2:end,:)=-D(:,:,:,2:2:end,:);
else
    D(:,:,:,1:2:end,:)=-D(:,:,:,1:2:end,:);
end
%% Quantify CBF (a anatomy featureless M0 is identified through quadratic interpolation)
Mult=6000*lambda*exp(PLD./T1b)./(2*alphaASL.*alphaASLInv*T1b.*(1-exp(-tau./T1b)));
%Mult=6000*lambda*exp(PLD./T1b)./(2*alpha.*alphaInv*T1b.*(1-exp(-tau./T1b)));

%     D=diff(ASL,1,4);%ASL(:,:,:,2:2:end,:)-ASL(:,:,:,1:2:end,:);
%     D(:,:,:,2:2:end,:)=-D(:,:,:,2:2:end,:);
%     lambda=0.9; %brain/blood partition coefficient;
%     PLD=1.5; % post label delay (s)
%     T1b=T1f*10^-3; % T1 of blood (s)
%     alpha=0.85; % tagging inversion efficiency
%     alphaInv=0.88; % efficiency correction due to background suppression
%     tau=1.5; % tagging time 1.5 s  
%     TR=4.4;
%     Mult=6000*lambda*exp(PLD./T1b)./(2*alpha.*alphaInv*T1b.*(1-exp(-tau./T1b)));

%Mult=6000*lambda*exp(PLD./T1b)./(2*alpha.*alphaInv*T1b.*(1-exp(-tau./T1b)));



M00=squeeze(M0b(:,:,:,1));
M00(M00==0)=NaN;
s=size(M00);
[X,Y,Z] = meshgrid(1:s(1),1:s(2),1:s(3));

%Set the basis functions
f1 = X.^2;
f2 = Y.^2;
f3 = Z.^2;
f4=Z.*Y;
f5=X.*Z;
f6=Y.*X;
f7=Y;
f8=X;
f9=Z;
f10=ones(size(X));
%Write as matrix equation
B = [f1(:),f2(:),f3(:), f4(:),f5(:),f6(:),f7(:),f8(:),f9(:),f10(:)];
y = M00(:);
%Solve for coefficients
B(isnan(y),:)=[];
y(isnan(y))=[];
indM=y>prctile(y,80);
indm=y<prctile(y,20);
y(logical(indM+indm))=[];
B(logical(indM+indm),:)=[];
c = B\y;

interpf=c(1)*f1+c(2)*f2+c(3)*f3+c(4)*f4+c(5)*f5+c(6)*f6+c(7)*f7+c(8)*f8+c(9)*f9+c(10)*f10;
M000=double(M00>0).*interpf;
M000(M000==0)=NaN;


m00=repmat(squeeze(M000),[1 1 1 size(D,4)]);
CBF=squeeze(D)./m00*Mult;
CBF(m00==0)=NaN;

%% Compute Average CBF map and Temporal SNR
CBF0=nanmean(CBF(:,:,:,baseline),4);
tSNR=CBF0./nanstd(CBF(:,:,:,baseline),[],4);

%% Estimate GM PVE and GM Mask from CBF map 
PVEGM=[];
GM=[];



PVEGM=CBF0/(prctile(CBF0(:),95)-prctile(CBF0(:),5));
GM=PVEGM>GMTH;

%% Extract Median GM Signal
CBFGM=[];

  
    for j=1:size(CBF,4)
        temp=squeeze(CBF(:,:,:,j));
        temp1=temp(:);
        CBFGM(j)=nanmedian(temp1(logical(GM(:))));
    end

DM0=mean(D,4).*double(M00>0);
DM0(DM0==0)=NaN;

ASL.CBF=CBF;
ASL.PERFUSION0=CBF0/Mult;
ASL.DM0=DM0;
ASL.CBF0=CBF0;
ASL.CBFGM=CBFGM;
ASL.PVEGM=PVEGM;
ASL.GM=GM;
ASL.tSNR=tSNR;
ASL.Parameters=Parameters;
%% Plot Data If Required
if FigOn==1
FigTitle=Parameters.FigTitle;
   FigH = figure('Position', get(0, 'Screensize'));
    for k=1:size(CBF0,3)

        subplot(5,5,k)
        imagesc(rot90(squeeze(ASL.DM0(:,:,k))))
        axis off
        axis equal
        caxis([0 prctile(ASL.DM0(:),95)])
    end
    set(gcf,'color','w');
    subplot(5,5,25)
    caxis([0 prctile(ASL.DM0(:),95)])
    colorbar('North')
    set(gca,'FontSize',12,'FontWeight','bold')
    axis off
    title('Baseline Perfusion Signal (A.U.)')
    subplot(5,5,23)
    axis off
    set(gca,'FontSize',12,'FontWeight','bold')
    subplot(5,5,3)
    title(FigTitle);
    axis off
    set(gca,'FontSize',12,'FontWeight','bold')



    FigH = figure('Position', get(0, 'Screensize'));
    for k=1:size(CBF0,3)

        subplot(5,5,k)
        imagesc(rot90(squeeze(100*ASL.PERFUSION0(:,:,k))))
        axis off
        axis equal
        caxis([0 prctile(100*ASL.PERFUSION0(:),95)])
    end
    set(gcf,'color','w');
    subplot(5,5,25)
    caxis([0 prctile(100*ASL.PERFUSION0(:),95)])
    colorbar('North')
    set(gca,'FontSize',12,'FontWeight','bold')
    axis off
    title('Baseline Perfusion Signal with M0 Normalization (%)')
    subplot(5,5,23)
    title(['Median GM Perfusion Signal=',num2str(round(100*nanmedian(100*ASL.PERFUSION0(GM)))/100),' %'])
    axis off
    set(gca,'FontSize',12,'FontWeight','bold')
    subplot(5,5,3)
    title(FigTitle);
    axis off
    set(gca,'FontSize',12,'FontWeight','bold')

    FigH = figure('Position', get(0, 'Screensize'));
    for k=1:size(CBF0,3)

        subplot(5,5,k)
        imagesc(rot90(squeeze(CBF0(:,:,k))))
        axis off
        axis equal
        caxis([0 prctile(CBF0(:),95)])
    end
    set(gcf,'color','w');
    subplot(5,5,25)
    caxis([0 prctile(CBF0(:),95)])
    colorbar('North')
    set(gca,'FontSize',12,'FontWeight','bold')
    axis off
    title('CBF0 (ml/100g/min)')
    subplot(5,5,23)
    title(['Median GM CBF0=',num2str(round(nanmedian(CBF0(GM)))),' ml/100g/min'])
    axis off
    set(gca,'FontSize',12,'FontWeight','bold')
    subplot(5,5,3)
    title(FigTitle);
    axis off
    set(gca,'FontSize',12,'FontWeight','bold')
   


    FigH = figure('Position', get(0, 'Screensize'));
    plot(1:length(ASL.CBFGM),ASL.CBFGM)
    xlabel('Samples')
    title([FigTitle,'; TR=',num2str(TR),' s'])
    set(gca,'FontSize',12,'FontWeight','bold')
    ylabel('Median GM CBF Timecourse (ml/100g/min)')
    set(gcf,'color','w');
    set(gca,'FontSize',12,'FontWeight','bold')



    FigH = figure('Position', get(0, 'Screensize'));
    for k=1:size(CBF0,3)
       
        subplot(5,5,k)
        imagesc(rot90(squeeze(tSNR(:,:,k))))
        axis off
        axis equal
        caxis([0 prctile(tSNR(:),95)])

    end
    set(gcf,'color','w');
    subplot(5,5,25)
    caxis([0 prctile(tSNR(:),95)])
    colorbar('North')
    set(gca,'FontSize',12,'FontWeight','bold')
    axis off
    title('Temporal SNR')
    subplot(5,5,23)
    temp1=tSNR(:);
    temp2=logical(GM(:));
    title(['CBF GM Median tSNR=',num2str(round(nanmedian(temp1(temp2)*100))/100)])
    axis off
    set(gca,'FontSize',12,'FontWeight','bold')
    subplot(5,5,3)
    title(FigTitle);
    axis off
    set(gca,'FontSize',12,'FontWeight','bold')


end


