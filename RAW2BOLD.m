% Written by Chiarelli Antonio Maria, University G. D'Annunzio of Chieti-Pescara
% Original Date: 12/03/2024
% Last Update: 05/03/2024




%The function process BOLD timecourses from BOLD-ASL fMRI.
% The volumes should be preprocessed including susceptibility distortions
% correction, motion correction and registration to M0
%Inputs:
% PathIn BOLD temporal volumes
% M0BrainPathIn: M0 brain masked map
% Parameters structure with the following field
%Parameters.baseline, baseline points where to evaluate average BOLD signal
%and temporal SNR
%Parameters.TR, Effective Repetition Time of dexi (seconds)
%Parameters.FigTitle= String with the title of the figures
%GM Mask of the Grey matter in the M0/ASL/BOLD space, set a 3D matrix with
%spatial dimensions equal to BOLD if GM not available
% FigOn, 1 plot figures, 0 do not plot figures

%Output:
%BOLD.BOLDraw, 1-sample moving average BOLD signal
%BOLD.BOLDavg, average BOLD signal
%BOLD.BOLDfra,fractional  BOLD signal changes
%BOLD.BOLDdet, fractional BOLD signal detrend
%BOLD.tSNR=tSNR; BOLD temporal SNR within the baseline period
%BOLD.BOLDGM=BOLDGM; average BOLD signal timecourses in the GM mask
%BOLD.GM=GM; GM mask



function [BOLD]=RAW2BOLD(PathIn,M0BrainPathIn,Parameters,GM,FigOn)

%% Import Images
TE2_hdr = spm_vol(PathIn);
TE2=spm_read_vols(TE2_hdr);
M0b_hdr = spm_vol(M0BrainPathIn);
M0b=spm_read_vols(M0b_hdr);

%% Import Parameters
TR=Parameters.TR; % Effective Repetition Time of dexi
FigTitle=Parameters.FigTitle;
baseline=Parameters.baseline;

%% Perform 1-sample moving average on BOLD (to remove the effect of tag-control alternation) and masking
BOLDraw=[];
for i=1:size(TE2,4)-1
   
    BOLDraw(:,:,:,i)=(TE2(:,:,:,i)+TE2(:,:,:,i+1))./2;
end
BOLDraw=BOLDraw.*repmat(M0b>0,[1 1 1 size(BOLDraw,4)]);
BOLDraw(BOLDraw==0)=NaN;
DIM=size(BOLDraw);

%% Compute Average BOLD map and Temporal SNR
BOLDavg=mean(BOLDraw(:,:,:,baseline),4);
tSNR=BOLDavg./nanstd(BOLDraw(:,:,:,baseline),[],4);

%% compute fractional BOLD and detrended BOLD
BOLDdet=NaN(DIM);
BOLDfra=(BOLDraw./repmat(BOLDavg,[1 1 1 DIM(4)])-1);
for ii=1:DIM(1)
    for jj=1:DIM(2)
        for kk=1:DIM(3)
            if ~isnan(tSNR(ii,jj,kk))
                signal=squeeze(BOLDfra(ii,jj,kk,:));
                BOLDdet(ii,jj,kk,:)= signal'-polyval(polyfit(1:length(signal),signal,1),1:length(signal));

            end
        end
    end
end

%% Extract Median GM Signal
BOLDGM=[];
for j=1:size(BOLDdet,4)
    temp=squeeze(BOLDdet(:,:,:,j));
    temp1=temp(:);
    BOLDGM(j)=nanmedian(temp1(logical(GM(:))));
end

%% save variables to the output structure

BOLD.BOLDraw=BOLDraw;
BOLD.BOLDavg=BOLDavg;
BOLD.BOLDdet=BOLDdet;
BOLD.BOLDfra=BOLDfra;
BOLD.tSNR=tSNR;
BOLD.BOLDGM=BOLDGM;
BOLD.GM=GM;
BOLD.Parameters=Parameters;

%% Plot Data If Required
if FigOn==1
    FigTitle=Parameters.FigTitle;
    FigH = figure('Position', get(0, 'Screensize'));
    for k=1:size(BOLDavg,3)

        subplot(5,5,k)
        imagesc(rot90(squeeze(BOLDavg(:,:,k))))
        axis off
        axis equal
        caxis([0 prctile(BOLDavg(:),95)])
    end
    set(gcf,'color','w');
    subplot(5,5,25)
    caxis([0 prctile(BOLDavg(:),95)])
    colorbar('North')
    set(gca,'FontSize',12,'FontWeight','bold')
    axis off
    title('Average BOLD Signal (A.U.)')
    subplot(5,5,3)
    title(FigTitle);
    axis off
    set(gca,'FontSize',12,'FontWeight','bold')



    FigH = figure('Position', get(0, 'Screensize'));
    plot(1:length(BOLD.BOLDGM),100*BOLD.BOLDGM)
    xlabel('Samples')
    title([FigTitle,'; TR=',num2str(TR),' s'])
    set(gca,'FontSize',12,'FontWeight','bold')
    ylabel('Median GM BOLD Timecourse (%)')
    set(gcf,'color','w');
    set(gca,'FontSize',12,'FontWeight','bold')



    FigH = figure('Position', get(0, 'Screensize'));
    for k=1:size(BOLDavg,3)

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
    title(['BOLD GM Median tSNR=',num2str(round(nanmedian(temp1(temp2)*100))/100)])
    axis off
    set(gca,'FontSize',12,'FontWeight','bold')
    subplot(5,5,3)
    title(FigTitle);
    axis off
    set(gca,'FontSize',12,'FontWeight','bold')


end


