
These set of functions use  concurrent BOLD-ASL fMRI recordings (raw volumes should be separated between BOLD and ASL and preprocessed with susceptibility distortion
and motion correction as well as registration to M0).
The functions process the data and, by inverting the Davis model of BOLD extended with the flow-diffusion
model of oxygen transport, estimate baseline OEF and CMRO2 by exploiting a vasodilatatory signal (single stimulus calibrated fMRI).
The vasodilatatory signal can be induced by hypercapnia, through CO2 inhalation or breath holding, or it can be related to an endogenous vascular signals (e.g., global BOLD signal or breathing signal in resting state data).
Please start from the Example_call_vasfMRI.m, using an exemplar BOLD-ASL dataset (generated with a breath holding task) to understand the different steps of the approach.
The different functions are described within each function header. 

Please refer to the following paper:
Chiarelli, A. M., Germuska, M., Chandler, H., Stickland, R., Patitucci, E., Biondetti, E., ... & Wise, R. G. (2022). A flow-diffusion model of oxygen transport for quantitative mapping of cerebral metabolic rate of oxygen (CMRO2) with single gas calibrated fMRI. Journal of Cerebral Blood Flow & Metabolism, 42(7), 1192-1209.


