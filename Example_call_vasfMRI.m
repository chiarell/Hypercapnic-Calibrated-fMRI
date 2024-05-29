 PathInCBF='/storage/shared/Code/CalibratedfMRI/SampleData/bh_asl_raw_mcf.nii.gz';
 PathInBOLD='/storage/shared/Code/CalibratedfMRI/SampleData/bh_bold_raw_mcf.nii.gz';
 M0BrainPathIn='/storage/shared/Code/CalibratedfMRI/SampleData/bh_M0_brain_ASL.nii.gz';
 Parameters.FigTitle = 'Sample Breath Hold';

Parameters.FigTitle = 'Sample Breath Hold';

Parameters.lambda=0.9; %brain/blood partition coefficient;
Parameters.PLD=1.5; % post label delay (s)
Parameters.T1b=1.54; % T1 of blood (s)
Parameters.alphaASL=0.85; % tagging inversion efficiency
Parameters.alphaASLInv=0.88; % efficiency correction due to background suppression
Parameters.tau=1.5; % tagging time 1.5 s
Parameters.TR=4.4; %  TR in s 
Parameters.TE=30; % TE in ms
Parameters.GMTH=0.5; % GM threashold
Parameters.baseline=1:121; % samples where to evaluate baseline values and temporal SNR
Parameters.HPorder=4;
Parameters.LPorder=8;
Parameters.TE=30;  %Time of echo of the BOLD sequence, units of ms
Parameters.TR=4.4; %Time of repetition, units of s
Parameters.THP=150; %high pass cut-off time, units of s
Parameters.TLP=10; %low pass cut-off time, units of s
Parameters.beta=1.3; % adimensional
Parameters.alpha=0.38; % Grubb exponent
Parameters.phi=1.32; % Oxygen binding macity of hemoglobin, ml of oxygen/g of Hb
Parameters.epsilon=0.0031; % Oxygen plasma solubility, ml/mmHg/dL
Parameters.h=2.8; % Hill coefficient, adimensional
Parameters.Hb=14;  %hemoglobin in blood g/dL
Parameters.OEFmin=0.01; % minimum OEF allowed
Parameters.OEFmax=1; % maximum OEF allowed
Parameters.OEFstep=10^-3; %OEF step for explicit search during the inversion (smaller steps slow the procedure, coarser steps make the estimate coarser)
Parameters.Mmin=0.01; % minimum M allowed
Parameters.Mmax=1; % maximum M allowed
Parameters.PmO2=0; % assumed average partial pressure at the mithocondria (mmHg)
%     %(units of mmHg, close to zero in Healthy subjects, i would use 10 in disease);
Parameters.Aro_k=8.8; % units of s^-1g^-beta dL^beta/(micromol/mmHg/ml/min) (suggested 10 s^-1g^-beta dL^beta/(micromol/mmHg/ml/min) at 3 T)
Parameters.n=10^5; %relative change in CBF over relative change in CMRO2
Parameters.shift=10; % shift for the regressor allowed, s
Parameters.InvMethod='Grid'; % inverse method, either Grid, Newton, or Search
%% Pre-process ASL to extract CBF

tag=0;
FigOn=1;
[ASL]=RAW2CBF(PathInCBF,M0BrainPathIn,tag,Parameters,FigOn);

%% Pre-process BOLD (fractional changes evaluation and detrending)
FigOn=1;
[BOLD]=RAW2BOLD(PathInBOLD,M0BrainPathIn,Parameters,ASL.GM,FigOn);

%% Filter CBF, BOLD and End-tidals
PETO2=randn(size(ASL.CBF,4),1)+100; % insert here endtidals in mmHg resampled at TR if available
PETCO2=randn(size(ASL.CBF,4),1)+36;
FigOn=1;
[fMRI_filt]=FilterCBF_BOLD_ET(ASL,BOLD,PETO2,PETCO2,Parameters,FigOn);

%% Estimate vas-fMRI Parameters
FigOn=1;
GMonly=0; %
RegressorType=2; %% normalized GM BOLD signal as a vasodilatatory signal
[outcome]=vasfMRI(fMRI_filt,GMonly,RegressorType,Parameters,FigOn);


%% Estimate Confidence in the hc-fMRI analysis UNDER CONSTRUCTION, WORK IN PROGRESS
n_iterations=10^4;
FigOn=1;
[outcome]=vasfMRIConfidence(outcome,fMRI_filt,ASL,n_iterations,FigOn);

