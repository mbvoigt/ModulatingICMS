%% Load TFR - AEES
%
% Subfunction for ModulatingICMS_main
%
% Loads TFRs for electrically stimulated trials
%
function [pow,plf] = Load_TFR_AEES(cfg)

% Initialize variables
pow = zeros(length(cfg.experiments),length(cfg.fi),length(cfg.ch2analyse),length(cfg.FTcfg.foi),1001);
plf = zeros(length(cfg.experiments),length(cfg.fi),length(cfg.ch2analyse),length(cfg.FTcfg.foi),1001);

for eID = 1:length(cfg.experiments)
    for f = 1:length(cfg.fi)
        
        fprintf('TEI-Analysis| %s | AEES | Stim. electrode %s | Loading TFR\n',cfg.experiments{eID},cfg.fi{f})
        mbv_log(sprintf('TEI-Analysis| %s | AEES | Stim. electrode %s | Loading TFR',cfg.experiments{eID},cfg.fi{f}),'file',cfg.logfile);
        % Load processed data file
        load(['C:\mbv\temp\tfr\TFR_AEES_' cfg.experiments{eID} '_' cfg.fi{f}]);
        
        % Get power spectrum
        pow(eID,f,1:length(cfg.ch2analyse),:,:) = TFR.powspctrm(cfg.ch2analyse,:,:);
        
        % Baseline subtraction
        for ch = 1:length(cfg.ch2analyse)
            for freq = 1:length(cfg.FTcfg.foi)
                baseline = squeeze(nanmean(pow(eID,f,ch,freq,cfg.BaselineWindow)));
                pow(eID,f,ch,freq,:) = 10*log10(pow(eID,f,ch,freq,:)/baseline);
            end
        end
        
        % Get phase-locking factor
        plf(eID,f,1:length(cfg.ch2analyse),:,:) = TFR.plfspctrm(cfg.ch2analyse,:,:);
        
        clear TFR
    end % repeat for every stimulated electrode
end % repeat for every experiment

% Remove channel 33 from experiments 1&2
% (Channel 33 was the current monitor resistor in experiments 1&2)
try
    if any(cfg.ch2analyse==33)
        pow([1,2],:,cfg.ch2analyse==33,:,:) = NaN;
    end
catch
    fprintf('Tried to remove channel 33 from experiments 1&2: Unsuccessful!\n')
end
