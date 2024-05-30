
These set of functions use  concurrent BOLD-ASL fMRI recordings (raw volumes should be separated between BOLD and ASL and preprocessed with susceptibility distortion
and motion correction as well as registration to M0).
The functions process the data and, by inverting the Davis model of BOLD extended with a flow-diffusion
model of oxygen transport, estimate baseline OEF and CMRO2 by exploiting a hypercapnic isometabolic modulation of blood flow and volume (single stimulus calibrated fMRI based on CO2 gas inhalation or breath holding).
The approach should work with any kind of flow modulation beyond hypercapnia with limitation generated just by the decreased SNR. The method, for example, can exploit endogenous modulations in brain hemodynamics (e.g., global BOLD signal or breathing pattern at rest).
Please start from the Example_call_vasfMRI.m, using an exemplar BOLD-ASL dataset (generated with a breath holding task) to understand the different steps of the approach.
The functions are described within each function header. 

Please refer to the following paper:
Chiarelli, A. M., Germuska, M., Chandler, H., Stickland, R., Patitucci, E., Biondetti, E., ... & Wise, R. G. (2022). A flow-diffusion model of oxygen transport for quantitative mapping of cerebral metabolic rate of oxygen (CMRO2) with single gas calibrated fMRI. Journal of Cerebral Blood Flow & Metabolism, 42(7), 1192-1209.


