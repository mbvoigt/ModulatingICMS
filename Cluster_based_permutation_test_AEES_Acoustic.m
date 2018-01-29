%% CBPT - Acoustic vs. AEES
%
% Subfunction for ModulatingICMS_main
%
% Performs the cluster-based permutation tests between the acoustically and 
% electrically stimulated condition 
%
% Plots figure 4C
%
function Cluster_based_permutation_test_AEES_Acoustic(cfg,pow_click,pow_AEES)

% Cluster-based permutation test: Configuration
cfg.CBPTcfg = [];
cfg.CBPTcfg.channel = 1:6;
cfg.CBPTcfg.frequency = 'all';
cfg.CBPTcfg.method = 'montecarlo';
cfg.CBPTcfg.statistic = 'ft_statfun_depsamplesT';
cfg.CBPTcfg.correctm = 'cluster';
cfg.CBPTcfg.clusteralpha = 0.05;
cfg.CBPTcfg.clusterstatistic = 'maxsum';
cfg.CBPTcfg.numrandomization = 1000;
cfg.CBPTcfg.design = [1 2 3 4 5 6 7 8 9 1 2 3 4 5 6 7 8 9;1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2];
cfg.CBPTcfg.uvar = 1;
cfg.CBPTcfg.ivar = 2;
cfg.CBPTcfg.neighbours = [];
cfg.CBPTcfg.tail = 0;
cfg.CBPTcfg.clustertail = 0;
cfg.CBPTcfg.alpha = 0.025;
cfg.CBPTcfg.minnbchan = 0;

% Cluster-based permutation test: Acoustic vs. AEES
% Shank 1
for stimL = 1:6
temp = squeeze(mean(pow_AEES(:,cfg.layerassign{stimL},:,:,:),2));
% Perform the test:
stat = Cluster_based_permutation_test(cfg,squeeze(nanmean(pow_click(:,1:16,:,:),2)),squeeze(nanmean(temp(:,1:16,:,:),2)));

% ---------
% FIGURE 4C - CBPT: Acoustiv vs. AEES
% ---------
fH = figure();
pcolor(cfg.FTcfg.toi,cfg.FTcfg.foi,squeeze(stat.stat(1,:,:)))
shading interp
hold on
contour(cfg.FTcfg.toi,cfg.FTcfg.foi,squeeze(stat.mask(1,:,:)),1,'Color','k','LineWidth',1.5)
caxis([-10 10])
colormap(jet)
print(fH,'-dpng','-r1800',['C:\mbv\temp\tfr\plots\CBPT_click_AEES_shank1_stimlayer' num2str(stimL)])
close(fH)
% ---------
end % repeat for each cortical layer stimulated

% Shank 2
for stimL = 1:6
temp = squeeze(mean(pow_AEES(:,cfg.layerassign{stimL},:,:,:),2));
% Perform the test:
stat = Cluster_based_permutation_test(cfg,squeeze(nanmean(pow_click(:,17:32,:,:),2)),squeeze(nanmean(temp(:,17:32,:,:),2)));

fH = figure();
pcolor(cfg.FTcfg.toi,cfg.FTcfg.foi,squeeze(stat.stat(1,:,:)))
shading interp
hold on
contour(cfg.FTcfg.toi,cfg.FTcfg.foi,squeeze(stat.mask(1,:,:)),1,'Color','k','LineWidth',1.5)
caxis([-10 10])
colormap(jet)
print(fH,'-dpng','-r1800',['C:\mbv\temp\tfr\plots\CBPT_click_AEES_shank2_stimlayer' num2str(stimL)])
close(fH)
end  % repeat for each cortical layer stimulated

