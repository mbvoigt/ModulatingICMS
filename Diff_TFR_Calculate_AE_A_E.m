%% Differential Analysis - Calculate AE - (A+E) TFRs
%
% Subfunction for ModulatingICMS_main
%
% Calculates the subtraction between combined (AE) and the sum of 
% acoustically (A) and electrically (E) stimulated TFRs:
% AE - (A+E)
%
%
function [res] = Diff_TFR_Calculate_AE_A_E(cfg,pow_AE,addA_E)

% Initialize variable
res = zeros(length(cfg.experiments),length(cfg.fi),length(cfg.ch2analyse),length(cfg.FTcfg.foi),1001,'single');

for eID = 1:length(cfg.experiments)
    for stimLayer = 1:length(cfg.fi)
        for ch = 1:length(cfg.ch2analyse)
            fprintf('%s | %02d | %02d\n',cfg.experiments{eID},stimLayer,ch)
            for freq = 1:length(cfg.FTcfg.foi)
                    res(eID,stimLayer,ch,freq,:) = squeeze(pow_AE(eID,stimLayer,cfg.ch2analyse(ch),freq,:))-squeeze(addA_E(eID,stimLayer,ch,freq,:));
            end % repeat for every frequency
        end % repeat for every channel
    end % repeat fpr every stimulated electrode
end % repeat for every experiment
