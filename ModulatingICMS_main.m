%% Intracortical microstimulation naturally modulates induced responses
%
% The present script, inlcuding all accompanying sub-scripts, was used to
% analyze the data and prepare the figures for the manuscript entitled
% "Intracortical microstimulation naturally modulates induced responses" by
% Mathias B. Voigt, Prasandhya A. Yusuf and Andrej Kral.
%
% Processing is arranged to mostly follow the occurence of the figures in 
% the manuscript. This is neither the order in which the script was run 
% from top to bottom, nor how it was written chronologically.
%
% Abbreviations:
%   AEES - All Electrodes Electrically Stimulated
%   CBPT - Cluster-based permutation test
%   TEI - Tone-Electric Interference
%   LFP - Local field potential
%   MUA - Multi unit activity
%
%
% External dependencies: 
%   - Matlab                                         (Version used: R2016b)
%   - Signal Processing Toolbox                      (Verison used: 7.3)
%   - Statistics and Machine Learning Toolboxes      (Version used: 11.0)
%   - fieldtrip                                      (Version used: r7276)
%
%
% This work is licensed under the MIT license.
% See the accompanying LICENSE.txt for detail.
%
% Copyright (c) 2018 Mathias Voigt <voigt.mathias@mh-hannover.de>
% Institute of AudioNeuroTechnology, Hannover Medical School
%
%
function ModulatingICMS_main

% Cleaning the workspace
clear variables
close all
clc

% Logging of script progress
cfg.logfile = ['C:\mbv\log\Publication_TEI_TFR_' datestr(now,'yy-mm-dd_HH_MM_SS') '.txt'];
mbv_log('Preparing TEI_TFR processing','file',cfg.logfile);
mbv_log('Preparing TEI_TFR processing','file',cfg.logfile);

% IDs of experiments to be analyzed
cfg.experiments = {'GP14F16CL','GP16F16CL','GP23F16CL','GP25J16CL','GP08K16CL','GP15K16CL','GP17K16CL','GP22K16CL','GP24K16CL'};
mbv_log(sprintf('Number of experiments to analyze: %2d',length(cfg.experiments)),'file',cfg.logfile);
mbv_log(['Experiments to analyze: ' strjoin(cfg.experiments,', ')],'file',cfg.logfile);
%
% Removed experiment with ID 'GP13I16CL', due to electrode placement in a 
% secondary auditory field!
%

% Acoustic stimulation was done with intensities of 120:-5:40 dB
% attenuation and coded as 01 to 17 (17 = 40 dB_Att)
% The following values are taken from the experimental protocols containing
% ABR threshold values to correspond to the data file containing the +40 dB
% intensity stimulation data for each experiment:
cfg.clickfile = {'12','11','10','13','11','11','11','12','14'};
%
% Removed the intensity corresponding to GP13I16CL ('11')
%

% Number of recorded channels 
cfg.noChannels = 48;
mbv_log(['Number of channels: ' num2str(cfg.noChannels)],'file',cfg.logfile);

%
% Channel mapping from the Neuronexus double-shank electrode array contacts
% to the Alpha Omega recording channels
%
cfg.channelidx=[9,8,10,7,11,6,12,5,13,4,14,3,15,2,16,1,...  % Shank 1
    25,24,26,23,27,22,28,21,29,20,30,19,31,18,32,...        % Shank 2
    17,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48];    % ECoG array
mbv_log(['Channel order used: ' sprintf('%02d, ',cfg.channelidx)],'file',cfg.logfile);

% Order of electrically stimulated channels
cfg.fi = {'08','07','09','06','10','05','11','04','12','03','13','02','14','01','15','00'};
mbv_log(['Stimulation channel order: ' strjoin(cfg.fi,', ')],'file',cfg.logfile);

% Stimulus delays (delta_t) analyzed
cfg.delay = {'5','15','25'};
mbv_log(['Delays [in ms] to analyze: ' strjoin(cfg.delay,', ')],'file',cfg.logfile);

% Time windows for linear interpolation to blank the electrical stimulation
% artefact, 'blf' blanking start time [ms], 'blt' blanking stop time [ms]
cfg.blf = [505;515;525];
cfg.blt = [508;518;528];
mbv_log(['Blanking interval starts: ' sprintf('%02d, ',cfg.blf)],'file',cfg.logfile);
mbv_log(['Blankin interval stops: ' sprintf('%02d, ',cfg.blt)],'file',cfg.logfile);

% Time window for epoch cutting of continuous recordings
% Here: 500 ms before and after stimulation trigger timestamp
cfg.sweepBefore = 500*22; 
cfg.sweepAfter = 500*22;
mbv_log('Analyzing 500 ms before and after trigger','file',cfg.logfile);

% Recalculate TFRs or use pre-existing ones
% [1 = new, 0 = use existing ones]
cfg.recalculate = 1;
if cfg.recalculate
    mbv_log('Recalculating TFRs and discarding existing ones!','file',cfg.logfile);
end

% Plot data figures or skip plotting 
% [1 = create figures, 0 = skip]
cfg.createPlots = 1;
if cfg.createPlots
    mbv_log('Creating data plots and discarding existing ones!','file',cfg.logfile);
end

% Fieldtrip configuration for Time-Frequency-Representation calculations
cfg.FTcfg.method = 'wavelet';
cfg.FTcfg.width = 6;
cfg.FTcfg.output = 'fourier';
cfg.FTcfg.channel = 1:48;
cfg.FTcfg.keeptrials = 'yes';
cfg.FTcfg.foi = 7:2:95;
cfg.FTcfg.toi = -0.5:0.001:0.5;

%-------------------
% Default parameters
cfg.ch2analyse = 1:cfg.noChannels;

% Time Windows [ms = samples@1kHz sampling rate]
cfg.EvokedResponse = 500:600;
cfg.InducedResponse = 600:800;
cfg.BaselineWindow = 1:351;

% Experiment assignment to the groups "Super-additive" and "Neutral"
% (Added for convenience for later analyses)
cfg.SupraAddExps = logical([1,1,0,0,0,0,0,1,1]);
cfg.NonRespoExps = ~cfg.SupraAddExps;

% Defininition of frequency bands
%
% alpha: < 15 Hz
% beta: 16 Hz - 30 Hz
% low gamma: 31 - 60 Hz
% high gamma: > 61 Hz 
%
cfg.alpha = 1:find(cfg.FTcfg.foi>15);
cfg.beta = find(cfg.FTcfg.foi>=16,1,'first'):find(cfg.FTcfg.foi>=30,1,'first');
cfg.logamma = find(cfg.FTcfg.foi>30,1,'first'):find(cfg.FTcfg.foi>=60,1,'first');
cfg.higamma = find(cfg.FTcfg.foi>60,1,'first'):length(cfg.FTcfg.foi);

% Electrode assignment to cortical layer according to Wallace&Palmer (2008)
cfg.layerassign = {1:2,3:4,5:6,7:8,9:12,13:16};
% Definition of supragranular layers (L2 & L3)
cfg.supgran = 3:6;

% Color definitions used for plotting
cfg.col_acoustic = [128/255 224/255 115/255];
cfg.col_acousticdark = [26/255 82/255 18/255];
cfg.col_electric = [255/255 140/255 0/255];
cfg.col_electricdark = [160/255 88/255 0/255]; 

%%
%
% -----------------------------------------------------------------------
% Local field potential (incl. ECoG) and multi unit activity calculations
% -----------------------------------------------------------------------
%
% Calculate Acoustic LFP response 
LFP_click = Calculate_LFP_Acoustic(cfg);

% Plot acoustically stimulated trials - local field potential 
Plot_LFP_Acoustic(cfg,LFP_click);
% Figures 2B, 3H, 3I

% Single-Trial-Phase analysis
Calculate_LFP_Acoustic_PhaseInfluence(cfg,LFP_click)
clear LFP_click
% Figures 7A, 7B, 7D, 7F

% Calculate Acoustic MUA response
MUA_click = Calculate_MUA_Acoustic(cfg);

% Plot acoustically stimulated trials - multi unit activity
Plot_MUA_Acoustic(cfg,MUA_click)
clear MUA_click
% Figure 2A

% Calculate AEES LFP response 
% (All Electrodes Electrically Stimulated = Electric only condition)
LFP_AEES = Calculate_LFP_AEES(cfg);

% Plot acoustically stimulated trials - local field potential
% (including model for surface amplitude when varying ICMS stimulation
% depth)
Plot_LFP_AEES(cfg,LFP_AEES);
clear LFP_AEES
% Figures 4G, 4H, 4I, 4J, 4K

% Single-Trial-Phase analysis
Calculate_LFP_AEES_PhaseInfluence(cfg)
% Figures 7C, 7E, 7G

%
% -----------------------------------------------------------------------
% Calculate Time-Frequency-Responses
% -----------------------------------------------------------------------
%
if cfg.recalculate
    Calculate_TFR_Acoustic(cfg)
    Calculate_TFR_AEES(cfg)
    Calculate_TFR_TEI(cfg)
end

%
% -----------------------------------------------------------------------
% Standard color bars
% -----------------------------------------------------------------------
%
% TFR power spectrum
fH = figure();
colormap(jet)
caxis([-5 5])
colorbar()
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\Colorbar_TFR_pow')
close(fH)
% TFR phase-locking factor
fH = figure();
colormap(jet)
caxis([0 1])
colorbar()
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\Colorbar_TFR_plf')
close(fH)


%
% -----------------------------------------------------------------------
% Acoustically stimulated TFRs
% -----------------------------------------------------------------------
%
% Load acoustic TFR data
[pow_click,plf_click] = Load_TFR_Acoustic(cfg);
% Plot acoustic TFR data
Plot_TFR_Acoustic(cfg,pow_click)
% Figure 2C
Plot_TFR_Acoustic(cfg,plf_click,'datatype','plf')
% Figure 3A

% Analyze TFR acoustic data
Analyze_TFR_Acoustic(cfg,pow_click,plf_click);
% Figures 3B, 3C, 3D, 3E, 3F, 3G
clear plf_click

%
% -----------------------------------------------------------------------
% Electrically stimulated TFRs
% -----------------------------------------------------------------------
%
% Load electric TFR data
[pow_AEES,plf_AEES] = Load_TFR_AEES(cfg);
% Plot electric TFR data
Plot_TFR_AEES(cfg,pow_AEES)
% Figure 4A
Plot_TFR_AEES(cfg,plf_AEES,'datatype','plf')
% Figure 4B
clear plf_AEES

% Cluster-based permutation test: Acoustic vs. AEES
Cluster_based_permutation_test_AEES_Acoustic(cfg,pow_click,pow_AEES)
% Figure 4C

% Analyze TFR AEES data
Analyse_TFR_AEES(cfg,pow_AEES);
% Figures 4D, 4E, 4F

%
% -----------------------------------------------------------------------
% Acoustically _AND_ electrically stimulated TFRs
% -----------------------------------------------------------------------
%

% Save baseline corrected TEI-TFR data to file
Load_TFR_TEI(cfg,'save2file',1);

% Seperately process data for each 16 electrode blocks:
% 1:16 = shank 1
% 17:32 = shank 2
% 33:48 = ECoG

% -------
% Shank 1
% -------
cfg.ch2analyse = 1:16;

% If variables for single stimulus modalities are not in workspace yet
% uncomment this:
% [pow_click,~] = Load_TFR_Acoustic(cfg);
% [pow_AEES,~] = Load_TFR_AEES(cfg);

% Differential Analysis - A+E
[~,addA_E] = Diff_TFR_Calculate_A_E(cfg,pow_click,pow_AEES);
save('C:\mbv\temp\tfr\addA_E_shank1','addA_E')
clear pow_AEES pow_click

% Load TEI data to workspace, calculate difference AE - (A+E) and save to
% file again

% delta_t = 5 ms
load('C:\mbv\temp\TEI_analyses_TEI_5ms_pow','pow_5ms')
res_5ms = Diff_TFR_Calculate_AE_A_E(cfg,pow_5ms,addA_E); %#ok variable saved to file
save('C:\mbv\temp\tfr\res_5ms_shank1','res_5ms')
clear pow_5ms res_5ms

% delta_t = 15 ms
load('C:\mbv\temp\TEI_analyses_TEI_15ms_pow','pow_15ms') 
res_15ms = Diff_TFR_Calculate_AE_A_E(cfg,pow_15ms,addA_E); %#ok variable saved to file
save('C:\mbv\temp\tfr\res_15ms_shank1','res_15ms')
clear pow_15ms res_15ms

% delta_t = 25 ms
load('C:\mbv\temp\TEI_analyses_TEI_25ms_pow','pow_25ms')
res_25ms = Diff_TFR_Calculate_AE_A_E(cfg,pow_25ms,addA_E); %#ok variable saved to file
save('C:\mbv\temp\tfr\res_25ms_shank1','res_25ms')
clear pow_25ms res_25ms

% -------
% Shank 2
% -------
clear addA_E

% Load shank 2 of acoustic and electric data
cfg.ch2analyse = 17:32;
[pow_click,~] = Load_TFR_Acoustic(cfg);
[pow_AEES,~] = Load_TFR_AEES(cfg);

cfg.ch2analyse = 1:16; % since only channels 17:32 are loaded only the 
% first 16 channels have to be analysed further 

% Differential Analysis - A+E
[~,addA_E] = Diff_TFR_Calculate_A_E(cfg,pow_click,pow_AEES);
save('C:\mbv\temp\tfr\addA_E_shank2','addA_E')
clear pow_AEES pow_click

cfg.ch2analyse = 17:32; % make sure shank 2 is processed for combined 
% stimulation!

% Load TEI data to workspace, calculate difference AE - (A+E) and save to
% file again

% delta_t = 5 ms
load('C:\mbv\temp\TEI_analyses_TEI_5ms_pow','pow_5ms')
res_5ms = Diff_TFR_Calculate_AE_A_E(cfg,pow_5ms,addA_E); %#ok variable saved to file
save('C:\mbv\temp\tfr\res_5ms_shank2','res_5ms')
clear pow_5ms res_5ms

% delta_t = 15 ms
load('C:\mbv\temp\TEI_analyses_TEI_15ms_pow','pow_15ms')
res_15ms = Diff_TFR_Calculate_AE_A_E(cfg,pow_15ms,addA_E); %#ok variable saved to file
save('C:\mbv\temp\tfr\res_15ms_shank2','res_15ms')
clear pow_15ms res_15ms

% delta_t = 25 ms
load('C:\mbv\temp\TEI_analyses_TEI_25ms_pow','pow_25ms')
res_25ms = Diff_TFR_Calculate_AE_A_E(cfg,pow_25ms,addA_E); %#ok variable saved to file
save('C:\mbv\temp\tfr\res_25ms_shank2','res_25ms')
clear pow_25ms res_25ms

% -------
%  ECoG
% -------
clear addA_E

% Load ECoG electrodes, acoustic and electric data
cfg.ch2analyse = 33:48;
[pow_click,~] = Load_TFR_Acoustic(cfg);
[pow_AEES,~] = Load_TFR_AEES(cfg);

cfg.ch2analyse = 1:16; % since only channels 33:48 are loaded only the 
% first 16 channels have to be analysed further 

% Differential Analysis - A+E
[~,addA_E] = Diff_TFR_Calculate_A_E(cfg,pow_click,pow_AEES);
save('C:\mbv\temp\tfr\addA_E_ECoG','addA_E')
clear pow_AEES pow_click

cfg.ch2analyse = 33:48; % make sure ECoG data is processed!

% Load TEI data to workspace, calculate difference AE - (A+E) and save to
% file again

% delta_t = 5 ms
load('C:\mbv\temp\TEI_analyses_TEI_5ms_pow','pow_5ms')
res_5ms = Diff_TFR_Calculate_AE_A_E(cfg,pow_5ms,addA_E); %#ok variable saved to file
save('C:\mbv\temp\tfr\res_5ms_ECoG','res_5ms')
clear pow_5ms res_5ms

% delta_t = 15 ms
load('C:\mbv\temp\TEI_analyses_TEI_15ms_pow','pow_15ms')
res_15ms = Diff_TFR_Calculate_AE_A_E(cfg,pow_15ms,addA_E); %#ok variable saved to file
save('C:\mbv\temp\tfr\res_15ms_ECoG','res_15ms')
clear pow_15ms res_15ms

% delta_t = 25 ms
load('C:\mbv\temp\TEI_analyses_TEI_25ms_pow','pow_25ms')
res_25ms = Diff_TFR_Calculate_AE_A_E(cfg,pow_25ms,addA_E); %#ok variable saved to file
save('C:\mbv\temp\tfr\res_25ms_ECoG','res_25ms')
clear pow_25ms res_25ms

clear addA_E

% ----------------------------------------
% Perform analyses on combined stimulation
% ----------------------------------------

Analyze_TFR_TEI(cfg)
% Plots figures 5A, 5B, 5C, 5D, 5E, 5F, 6A, 6B, 6C, 6D, 6E, 6F

% Single-Trial-Phase analysis
Calculate_LFP_TEI_PhaseInfluence(cfg)
% Plots figures 7H, 7I, 7J


% ----------------------------------------
% ----------------------------------------
% ----------------------------------------
fprintf('TEI_TFR processing DONE!\n');
mbv_log('TEI_TFR processing DONE!','file',cfg.logfile);

end