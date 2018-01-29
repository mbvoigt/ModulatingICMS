%% Load TFR - Acoustic
%
% Subfunction for ModulatingICMS_main
%
% Loads TFRs for acoustically stimulated trials
%
function [pow,plf] = Load_TFR_Acoustic(cfg)

% Initialize variables
pow = zeros(length(cfg.experiments),length(cfg.ch2analyse),length(cfg.FTcfg.foi),1001);
plf = zeros(length(cfg.experiments),length(cfg.ch2analyse),length(cfg.FTcfg.foi),1001);

for eID = 1:length(cfg.experiments)
    
    fprintf('TEI-Analysis| %s | Acoustic | Loading TFR\n',cfg.experiments{eID})
    mbv_log(sprintf('TEI-Analysis| %s | Acoustic | Loading TFR',cfg.experiments{eID}),'file',cfg.logfile);
    % Load processed data file
    load(['C:\mbv\temp\tfr\TFR_' cfg.experiments{eID} '_click']);
    
    % Get power spectrum
    pow(eID,1:length(cfg.ch2analyse),:,:) = TFR.powspctrm(cfg.ch2analyse,:,:);
    
    % Baseline subtraction
    for ch = 1:length(cfg.ch2analyse)
        for freq = 1:length(cfg.FTcfg.foi)
            baseline = squeeze(nanmean(pow(eID,ch,freq,cfg.BaselineWindow)));
            pow(eID,ch,freq,:) = 10*log10(pow(eID,ch,freq,:)/baseline);
        end
    end
    
    % Get phase-locking factor
    plf(eID,1:length(cfg.ch2analyse),:,:) = TFR.plfspctrm(cfg.ch2analyse,:,:);
    
    clear TFR
end % repeat for every experiment

% Remove channel 33 from experiments 1&2
% (Channel 33 was the current monitor resistor in experiments 1&2)
try
    if any(cfg.ch2analyse==33)
        pow([1,2],cfg.ch2analyse==33,:,:) = NaN;
    end
catch
    fprintf('Tried to remove channel 33 from experiments 1&2: Unsuccessful!\n')
end
