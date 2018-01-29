%% Differential Analysis - Calculate A+E TFRs
%
% Subfunction for ModulatingICMS_main
%
% Calculates the addition and subtraction between acoustically and
% electrically stimulated TFRs = Acoustic -/+ AEES
%
function [difA_E, addA_E] = Diff_TFR_Calculate_A_E(cfg,pow_click,pow_AEES)

% Initialize variablse
difA_E = zeros(length(cfg.experiments),length(cfg.fi),length(cfg.ch2analyse),length(cfg.FTcfg.foi),1001);
addA_E = zeros(length(cfg.experiments),length(cfg.fi),length(cfg.ch2analyse),length(cfg.FTcfg.foi),1001);

for eID = 1:length(cfg.experiments)
    for stimLayer = 1:length(cfg.fi)
        for ch = 1:length(cfg.ch2analyse)
            fprintf('%s | %02d | %02d\n',cfg.experiments{eID},stimLayer,cfg.ch2analyse(ch))
            for freq = 1:length(cfg.FTcfg.foi)
                difA_E(eID,stimLayer,ch,freq,:) = squeeze(pow_click(eID,cfg.ch2analyse(ch),freq,:))-squeeze(pow_AEES(eID,stimLayer,cfg.ch2analyse(ch),freq,:));
                addA_E(eID,stimLayer,ch,freq,:) = squeeze(pow_click(eID,cfg.ch2analyse(ch),freq,:))+squeeze(pow_AEES(eID,stimLayer,cfg.ch2analyse(ch),freq,:));
            end % repeat for each frequency
        end % repeat for each channel
    end % repeat for each stimulated electrode
end % repeat for each experiment
