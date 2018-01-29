%% Analyse TFR - Acoustic
%
% Subfunction for ModulatingICMS_main
%
% Analyzed TFRs for acoustically stimulated trials & plots further figures
%
% Plots figures 3B, 3C, 3D, 3E, 3F, 3G
%
function Analyze_TFR_Acoustic(cfg,pow_click,plf_click)

%% Phase-locking factor analysis
% Shank 1
plftmp = squeeze(nanmean(nanmean(plf_click(:,1:16,:,:),2),3));

% ---------
% FIGURE 3B - Phase-locking factor collapse
% ---------
fH=figure();
plot(nanmean(plftmp),'k') % mean over experiments
hold on
plot(nanmean(plftmp)+(nanstd(plftmp)/sqrt(9)),'--k') % + SEM
plot(nanmean(plftmp)-(nanstd(plftmp)/sqrt(9)),'--k') % - SEM
set(gca,'YLim',[0 1],'XLim',[0 1000],'XTick',0:100:1000,'XTickLabel',-0.5:0.1:0.5)
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_click_plf_shank1_collapse')
close(fH)
% ---------

%% Evoked & induced amplitude determination

% Initialize variables 
evokedAmplitude = NaN(length(cfg.experiments),cfg.noChannels,length(cfg.FTcfg.foi));
inducedAmplitude = NaN(length(cfg.experiments),cfg.noChannels,length(cfg.FTcfg.foi));
% Determine amplitude
for eID = 1:length(cfg.experiments)
    for ch = 1: cfg.noChannels
        evokedAmplitude(eID,ch,:) = mean(pow_click(eID,ch,:,cfg.EvokedResponse),4);
        inducedAmplitude(eID,ch,:) = mean(pow_click(eID,ch,:,cfg.InducedResponse),4);
    end
end

% ---------
% FIGURE 3C - Evoked response as a function of cortical layer
% ---------
% Shank 1
tmp = [max(max(evokedAmplitude(:,cfg.layerassign{1},:),[],3),[],2),...
    max(max(evokedAmplitude(:,cfg.layerassign{2},:),[],3),[],2),...
    max(max(evokedAmplitude(:,cfg.layerassign{3},:),[],3),[],2),...
    max(max(evokedAmplitude(:,cfg.layerassign{4},:),[],3),[],2),...
    max(max(evokedAmplitude(:,cfg.layerassign{5},:),[],3),[],2),...
    max(max(evokedAmplitude(:,cfg.layerassign{6},:),[],3),[],2)];
fH = figure();
boxplot(tmp,'Colors','k');
hold on
plot(tmp','o','MarkerFaceColor',cfg.col_acoustic,'MarkerEdgeColor','none')
set(gca,'YLim',[0 25])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_click_TFR_evokedAmp_layers_shank1')
close(fH)
% ---------

% Statistics
[p,tbl,stats]=anova1(tmp); %#ok Results used for manuscript text but not for further calculations

% Shank 2
tmp = [max(max(evokedAmplitude(:,cfg.layerassign{1}+16,:),[],3),[],2),...
    max(max(evokedAmplitude(:,cfg.layerassign{2}+16,:),[],3),[],2),...
    max(max(evokedAmplitude(:,cfg.layerassign{3}+16,:),[],3),[],2),...
    max(max(evokedAmplitude(:,cfg.layerassign{4}+16,:),[],3),[],2),...
    max(max(evokedAmplitude(:,cfg.layerassign{5}+16,:),[],3),[],2),...
    max(max(evokedAmplitude(:,cfg.layerassign{6}+16,:),[],3),[],2)];
fH = figure();
boxplot(tmp,'Colors','k');
hold on
plot(tmp','o','MarkerFaceColor',cfg.col_acoustic,'MarkerEdgeColor','none')
set(gca,'YLim',[0 25])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_click_TFR_evokedAmp_layers_shank2')
close(fH)
% ---------

% Statistics
[p,tbl,stats]=anova1(tmp); %#ok Results used for manuscript text but not for further calculations


% ---------
% FIGURE 3E - Evoked response as a function of frequency band
% ---------
% Shank 1
tmp = [max(max(evokedAmplitude(:,1:16,cfg.alpha),[],3),[],2),...
    max(max(evokedAmplitude(:,1:16,cfg.beta),[],3),[],2),...
    max(max(evokedAmplitude(:,1:16,cfg.logamma),[],3),[],2),...
    max(max(evokedAmplitude(:,1:16,cfg.higamma),[],3),[],2)];
fH = figure();
boxplot(tmp,'Colors','k');
hold on
plot(tmp','o','MarkerFaceColor',cfg.col_acoustic,'MarkerEdgeColor','none')
set(gca,'YLim',[0 25])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_click_TFR_evokedAmp_freqbands_shank1')
close(fH)
% ---------

% Statistics
[p,tbl,stats]=anova1(tmp);  %#ok Results used for manuscript text but not for further calculations

% Shank 2
tmp = [max(max(evokedAmplitude(:,17:32,cfg.alpha),[],3),[],2),...
    max(max(evokedAmplitude(:,17:32,cfg.beta),[],3),[],2),...
    max(max(evokedAmplitude(:,17:32,cfg.logamma),[],3),[],2),...
    max(max(evokedAmplitude(:,17:32,cfg.higamma),[],3),[],2)];
fH = figure();
boxplot(tmp,'Colors','k');
hold on
plot(tmp','o','MarkerFaceColor',cfg.col_acoustic,'MarkerEdgeColor','none')
set(gca,'YLim',[0 25])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_click_TFR_evokedAmp_freqbands_shank2')
close(fH)
% ---------

% Statistics
[p,tbl,stats]=anova1(tmp);  %#ok Results used for manuscript text but not for further calculations

% ---------
% FIGURE 3G - ECoG - Evoked response as a function of frequency band
% ---------
tmp = [max(max(evokedAmplitude(:,33:48,cfg.alpha),[],3),[],2),...
    max(max(evokedAmplitude(:,33:48,cfg.beta),[],3),[],2),...
    max(max(evokedAmplitude(:,33:48,cfg.logamma),[],3),[],2),...
    max(max(evokedAmplitude(:,33:48,cfg.higamma),[],3),[],2)];
fH = figure();
boxplot(tmp,'Colors','k');
hold on
plot(tmp','o','MarkerFaceColor',cfg.col_acoustic,'MarkerEdgeColor','none')
set(gca,'YLim',[0 25])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_click_TFR_evokedAmp_freqbands_ECoG')
close(fH)
% ---------

% Statistics
[p,tbl,stats]=anova1(tmp); %#ok Results used for manuscript text but not for further calculations

% ---------
% FIGURE 3D - Induced response as a function of cortical layer
% ---------
% Shank 1
tmp = [max(max(inducedAmplitude(:,cfg.layerassign{1},:),[],3),[],2),...
    max(max(inducedAmplitude(:,cfg.layerassign{2},:),[],3),[],2),...
    max(max(inducedAmplitude(:,cfg.layerassign{3},:),[],3),[],2),...
    max(max(inducedAmplitude(:,cfg.layerassign{4},:),[],3),[],2),...
    max(max(inducedAmplitude(:,cfg.layerassign{5},:),[],3),[],2),...
    max(max(inducedAmplitude(:,cfg.layerassign{6},:),[],3),[],2)];
fH = figure();
boxplot(tmp,'Colors','k');
hold on
plot(tmp','o','MarkerFaceColor',cfg.col_acoustic,'MarkerEdgeColor','none')
set(gca,'YLim',[0 25])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_click_TFR_inducedAmp_layers_shank1')
close(fH)
% ---------

% Statistics
[p,tbl,stats]=anova1(tmp); %#ok Results used for manuscript text but not for further calculations

% Shank 2
tmp = [max(max(inducedAmplitude(:,cfg.layerassign{1}+16,:),[],3),[],2),...
    max(max(inducedAmplitude(:,cfg.layerassign{2}+16,:),[],3),[],2),...
    max(max(inducedAmplitude(:,cfg.layerassign{3}+16,:),[],3),[],2),...
    max(max(inducedAmplitude(:,cfg.layerassign{4}+16,:),[],3),[],2),...
    max(max(inducedAmplitude(:,cfg.layerassign{5}+16,:),[],3),[],2),...
    max(max(inducedAmplitude(:,cfg.layerassign{6}+16,:),[],3),[],2)];
fH = figure();
boxplot(tmp,'Colors','k');
hold on
plot(tmp','o','MarkerFaceColor',cfg.col_acoustic,'MarkerEdgeColor','none')
set(gca,'YLim',[0 25])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_click_TFR_inducedAmp_layers_shank2')
close(fH)
% ---------

% Statistics
[p,tbl,stats]=anova1(tmp); %#ok Results used for manuscript text but not for further calculations


% ---------
% FIGURE 3F - Induced response as a function of frequency band
% ---------
% Shank 1
tmp = [max(max(inducedAmplitude(:,1:16,cfg.alpha),[],3),[],2),...
    max(max(inducedAmplitude(:,1:16,cfg.beta),[],3),[],2),...
    max(max(inducedAmplitude(:,1:16,cfg.logamma),[],3),[],2),...
    max(max(inducedAmplitude(:,1:16,cfg.higamma),[],3),[],2)];
fH = figure();
boxplot(tmp,'Colors','k');
hold on
plot(tmp','o','MarkerFaceColor',cfg.col_acoustic,'MarkerEdgeColor','w')
set(gca,'YLim',[0 25])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_click_TFR_inducedAmp_freqbands_shank1')
close(fH)
% ---------

% Statistics
[p,tbl,stats]=anova1(tmp); %#ok Results used for manuscript text but not for further calculations

% Shank 2
tmp = [max(max(inducedAmplitude(:,17:32,cfg.alpha),[],3),[],2),...
    max(max(inducedAmplitude(:,17:32,cfg.beta),[],3),[],2),...
    max(max(inducedAmplitude(:,17:32,cfg.logamma),[],3),[],2),...
    max(max(inducedAmplitude(:,17:32,cfg.higamma),[],3),[],2)];
fH = figure();
boxplot(tmp,'Colors','k');
hold on
plot(tmp','o','MarkerFaceColor',cfg.col_acoustic,'MarkerEdgeColor','w')
set(gca,'YLim',[0 25])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_click_TFR_inducedAmp_freqbands_shank2')
close(fH)
% ---------

% Statistics
[p,tbl,stats]=anova1(tmp); %#ok Results used for manuscript text but not for further calculations

% ---------
% FIGURE 3G - ECoG - Induced response as a function of frequency band
% ---------
tmp = [max(max(inducedAmplitude(:,33:48,cfg.alpha),[],3),[],2),...
    max(max(inducedAmplitude(:,33:48,cfg.beta),[],3),[],2),...
    max(max(inducedAmplitude(:,33:48,cfg.logamma),[],3),[],2),...
    max(max(inducedAmplitude(:,33:48,cfg.higamma),[],3),[],2)];
fH = figure();
boxplot(tmp,'Colors','k');
hold on
plot(tmp','o','MarkerFaceColor',cfg.col_acoustic,'MarkerEdgeColor','w')
set(gca,'YLim',[0 25])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_click_TFR_inducedAmp_freqbands_ECoG')
close(fH)
% ---------

% Statistics
[p,tbl,stats]=anova1(tmp); %#ok Results used for manuscript text but not for further calculations



