%% Cluster based permutation test
%
% Subfunction for
% ModulatingICMS_main/Cluster_based_permutation_test_AEES_Acoustic
%
% Performs CBPT using fieldtrip
%
function stat = Cluster_based_permutation_test(cfg,pow_1,pow_2)

TFR1 = struct;
TFR2 = struct;

% Reformat data to match format expected by fieldtrip CBPT 
for eID = 1:size(pow_1,1)
    TFR1{eID}.powspctrm(1,:,:) = squeeze(pow_1(eID,:,:));
    TFR1{eID}.dimord = 'chan_freq_time';
    TFR1{eID}.label = {'1'};
    TFR1{eID}.freq = cfg.FTcfg.foi;
    TFR1{eID}.time = cfg.FTcfg.toi;
    TFR1{eID}.cumtapcnt = ones(30,length(cfg.FTcfg.foi));
    %
    TFR2{eID}.powspctrm(1,:,:) = squeeze(pow_2(eID,:,:));
    TFR2{eID}.dimord = 'chan_freq_time';
    TFR2{eID}.label = {'1'};
    TFR2{eID}.freq = cfg.FTcfg.foi;
    TFR2{eID}.time = cfg.FTcfg.toi;
    TFR2{eID}.cumtapcnt = ones(30,length(cfg.FTcfg.foi));
end

GAcfg = [];
GAcfg.channel = 1;
GAcfg.keepindividual = 'yes';
GA1 = ft_freqgrandaverage(GAcfg,TFR1{:});
GA2 = ft_freqgrandaverage(GAcfg,TFR2{:});

% TEST:
[stat] = ft_freqstatistics(cfg.CBPTcfg,GA1,GA2);

end
