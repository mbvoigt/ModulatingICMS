%% Analyse TFR - TEI - Combined acoustic _and_ electric stimulation
%
% Subfunction for ModulatingICMS_main
%
% Analyzes TFRs for combined acoustic _and_ electric stimulated trials, 
% analyzes the supra-additive effects & plots further figures for both
%
% Plots figures 5A, 5B, 5C, 5D, 5E, 5F, 6A, 6B, 6C, 6D, 6E, 6F
%
function Analyse_TFR_TEI(cfg)

%% Difference calculation example
%
% ---------
% Figure 5A
% ---------
%
% Pictures to show calculation of supra-additive effect of combined
% stimulation
%

exampleExperiment = 1;
% more or less arbitrarily chosen from one of the experiments classified as
% super-additive (see below!)

% Acoustic total power example
cfg.ch2analyse = 1:16; % make sure electrodes 1:16 are loaded
[pow_click,~] = Load_TFR_Acoustic(cfg);

fH = figure();
pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(pow_click(exampleExperiment,cfg.supgran,:,:),2)))
shading interp
caxis([-5 5])
colormap(jet)
print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\DiffExample_click')
close(fH)

clear pow_click

% Electric total power example
[pow_AEES,~] = Load_TFR_AEES(cfg);

fH = figure();
pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(pow_AEES(exampleExperiment,cfg.supgran,cfg.supgran,:,:),2),3)))
shading interp
caxis([-5 5])
colormap(jet)
print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\DiffExample_AEES')
close(fH)

clear pow_AEES

% Combined total power example
load('C:\mbv\temp\TEI_analyses_TEI_5ms_pow','pow_5ms')

fH = figure();
pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(pow_5ms(exampleExperiment,cfg.supgran,cfg.supgran,:,:),2),3))) %#ok variable loaded from file
shading interp
caxis([-5 5])
colormap(jet)
print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\DiffExample_5msPower')
close(fH)

clear pow_5ms

% Difference example
load('C:\mbv\temp\tfr\res_5ms_shank1','res_5ms')

fH = figure();
pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(res_5ms(exampleExperiment,cfg.supgran,cfg.supgran,:,:),2),3))) %#ok variable loaded from file
shading interp
caxis([-5 5])
colormap(jet)
print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\DiffExample_5msDifference')
close(fH)

clear res_5ms exampleExperiment

%% Phase-locking factor example

% Initialize variable
plf_5ms = zeros(length(cfg.experiments),length(cfg.fi),48,45,1001);

% Load TEI PLF data
for eID = 1:length(cfg.experiments)
    for f = 1:length(cfg.fi)
        fprintf('%2d | %2d\n',eID,f)
        load(['C:\mbv\temp\tfr\TFR_' cfg.experiments{eID} '_5_' cfg.fi{f}]);
        % Get phase-locking factor
        plf_5ms(eID,f,:,:,:) = TFR.plfspctrm(:,:,:);
        clear TFR ST_LGP ST_LGPE
        % Remove stimulation channel
        plf_5ms(eID,f,f,:,:) = NaN;
    end
end
clear eID f

% ---------
% Figure 5B
% ---------
%
fH = figure();
pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(nanmean(plf_5ms(:,:,1:16,:,:),1),2),3))) % PLF of shank 1
shading interp
caxis([0 1])
colormap(jet)
print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\TEI_5ms_plf_shank1')
close(fH)

% ---------
% Figure 5C
% ---------
%
plftmp = squeeze(nanmean(nanmean(nanmean(plf_5ms(:,:,1:16,:,:),2),3),4)); % PLF of shank 1

fH=figure();
plot(nanmean(plftmp),'k') % mean
hold on
plot(nanmean(plftmp)+(nanstd(plftmp)/sqrt(9)),'--k') % +SEM
plot(nanmean(plftmp)-(nanstd(plftmp)/sqrt(9)),'--k') % -SEM
set(gca,'YLim',[0 1],'XLim',[0 1000],'XTick',0:100:1000,'XTickLabel',-0.5:0.1:0.5)
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_5ms_plf_shank1_collapse')
close(fH)

clear plf_5ms

%% Total power analysis

% Colormap for single lines
cm = flipud([20/255 11/255 0/255;...
    100/255 55/255 0/255;...
    160/255 88/255 0/255;...
    255/255 140/255 0/255;...
    255/255 160/255 45/255;...
    255/255 205/255 145/255]);

% -----
%  5 ms
% -----

load('C:\mbv\temp\TEI_analyses_TEI_5ms_pow','pow_5ms')

% Initialize variables
evAmp_5ms_1 = zeros(length(cfg.experiments),length(cfg,fi),45);
evAmp_5ms_2 = zeros(length(cfg.experiments),length(cfg,fi),45);
evAmp_5ms_3 = zeros(length(cfg.experiments),length(cfg,fi),45);
inAmp_5ms_1 = zeros(length(cfg.experiments),length(cfg,fi),45);
inAmp_5ms_2 = zeros(length(cfg.experiments),length(cfg,fi),45);
inAmp_5ms_3 = zeros(length(cfg.experiments),length(cfg,fi),45);

% Get evoked amplitude
for eID = 1:length(cfg.experiments)
    for f = 1:16
        evAmp_5ms_1(eID,f,:) = nanmean(max(pow_5ms(eID,f,cfg.supgran,:,cfg.EvokedResponse),[],5),3);
        evAmp_5ms_2(eID,f,:) = nanmean(max(pow_5ms(eID,f,cfg.supgran+16,:,cfg.EvokedResponse),[],5),3);
        evAmp_5ms_3(eID,f,:) = nanmean(max(pow_5ms(eID,f,33:48,:,500:600),[],5),3);
    end
end

% Get induced amplitude
for eID = 1:length(cfg.experiments)
    for f = 1:16
        inAmp_5ms_1(eID,f,:) = nanmean(max(pow_5ms(eID,f,cfg.supgran,:,cfg.InducedResponse),[],5),3);
        inAmp_5ms_2(eID,f,:) = nanmean(max(pow_5ms(eID,f,cfg.supgran+16,:,cfg.InducedResponse),[],5),3);
        inAmp_5ms_3(eID,f,:) = nanmean(max(pow_5ms(eID,f,33:48,:,cfg.InducedResponse),[],5),3);
    end
end

clear pow_5ms

% -----
% 15 ms
% -----

load('C:\mbv\temp\TEI_analyses_TEI_15ms_pow','pow_15ms')

% Initialize variables
evAmp_15ms_1 = zeros(length(cfg.experiments),length(cfg,fi),45);
evAmp_15ms_2 = zeros(length(cfg.experiments),length(cfg,fi),45);
evAmp_15ms_3 = zeros(length(cfg.experiments),length(cfg,fi),45);
inAmp_15ms_1 = zeros(length(cfg.experiments),length(cfg,fi),45);
inAmp_15ms_2 = zeros(length(cfg.experiments),length(cfg,fi),45);
inAmp_15ms_3 = zeros(length(cfg.experiments),length(cfg,fi),45);

% Get evoked amplitude
for eID = 1:length(cfg.experiments)
    for f = 1:16
        evAmp_15ms_1(eID,f,:) = nanmean(max(pow_15ms(eID,f,cfg.supgran,:,cfg.EvokedResponse),[],5),3); %#ok variable loaded from file
        evAmp_15ms_2(eID,f,:) = nanmean(max(pow_15ms(eID,f,cfg.supgran+16,:,cfg.EvokedResponse),[],5),3);
        evAmp_15ms_3(eID,f,:) = nanmean(max(pow_15ms(eID,f,33:48,:,cfg.EvokedResponse),[],5),3);
    end
end

% Get induced amplitude
for eID = 1:length(cfg.experiments)
    for f = 1:16
        inAmp_15ms_1(eID,f,:) = nanmean(max(pow_15ms(eID,f,cfg.supgran,:,cfg.InducedResponse),[],5),3);
        inAmp_15ms_2(eID,f,:) = nanmean(max(pow_15ms(eID,f,cfg.supgran+16,:,cfg.InducedResponse),[],5),3);
        inAmp_15ms_3(eID,f,:) = nanmean(max(pow_15ms(eID,f,33:48,:,cfg.InducedResponse),[],5),3);
    end
end

clear pow_15ms

% -----
% 25 ms
% -----

load('C:\mbv\temp\TEI_analyses_TEI_25ms_pow','pow_25ms')

% Initialize variables
evAmp_25ms_1 = zeros(length(cfg.experiments),length(cfg,fi),45);
evAmp_25ms_2 = zeros(length(cfg.experiments),length(cfg,fi),45);
evAmp_25ms_3 = zeros(length(cfg.experiments),length(cfg,fi),45);
inAmp_25ms_1 = zeros(length(cfg.experiments),length(cfg,fi),45);
inAmp_25ms_2 = zeros(length(cfg.experiments),length(cfg,fi),45);
inAmp_25ms_3 = zeros(length(cfg.experiments),length(cfg,fi),45);

% Get evoked amplitude
for eID = 1:length(cfg.experiments)
    for f = 1:16
        evAmp_25ms_1(eID,f,:) = nanmean(max(pow_25ms(eID,f,cfg.supgran,:,cfg.EvokedResponse),[],5),3); %#ok variable loaded from file
        evAmp_25ms_2(eID,f,:) = nanmean(max(pow_25ms(eID,f,cfg.supgran+16,:,cfg.EvokedResponse),[],5),3);
        evAmp_25ms_3(eID,f,:) = nanmean(max(pow_25ms(eID,f,33:48,:,cfg.EvokedResponse),[],5),3);
    end
end

% Get induced amplitude
for eID = 1:length(cfg.experiments)
    for f = 1:16
        inAmp_25ms_1(eID,f,:) = nanmean(max(pow_25ms(eID,f,cfg.supgran,:,cfg.InducedResponse),[],5),3);
        inAmp_25ms_2(eID,f,:) = nanmean(max(pow_25ms(eID,f,cfg.supgran+16,:,cfg.InducedResponse),[],5),3);
        inAmp_25ms_3(eID,f,:) = nanmean(max(pow_25ms(eID,f,33:48,:,cfg.InducedResponse),[],5),3);
    end
end

clear pow_25ms eID f

%%
% ---------
% Figure 5D
% ---------
%
% Shank 1
%
% Mean
tmp = [nanmean(nanmean(nanmean(evAmp_5ms_1(:,cfg.layerassign{1},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_1(:,cfg.layerassign{1},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_1(:,cfg.layerassign{1},:),2),3));...
    nanmean(nanmean(nanmean(evAmp_5ms_1(:,cfg.layerassign{2},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_1(:,cfg.layerassign{2},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_1(:,cfg.layerassign{2},:),2),3));...
    nanmean(nanmean(nanmean(evAmp_5ms_1(:,cfg.layerassign{3},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_1(:,cfg.layerassign{3},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_1(:,cfg.layerassign{3},:),2),3));...
    nanmean(nanmean(nanmean(evAmp_5ms_1(:,cfg.layerassign{4},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_1(:,cfg.layerassign{4},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_1(:,cfg.layerassign{4},:),2),3));...
    nanmean(nanmean(nanmean(evAmp_5ms_1(:,cfg.layerassign{5},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_1(:,cfg.layerassign{5},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_1(:,cfg.layerassign{5},:),2),3));...
    nanmean(nanmean(nanmean(evAmp_5ms_1(:,cfg.layerassign{6},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_1(:,cfg.layerassign{6},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_1(:,cfg.layerassign{6},:),2),3))];

% Standard deviation
tmps = [nanstd(nanmean(nanmean(evAmp_5ms_1(:,cfg.layerassign{1},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_1(:,cfg.layerassign{1},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_1(:,cfg.layerassign{1},:),2),3));...
    nanstd(nanmean(nanmean(evAmp_5ms_1(:,cfg.layerassign{2},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_1(:,cfg.layerassign{2},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_1(:,cfg.layerassign{2},:),2),3));...
    nanstd(nanmean(nanmean(evAmp_5ms_1(:,cfg.layerassign{3},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_1(:,cfg.layerassign{3},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_1(:,cfg.layerassign{3},:),2),3));...
    nanstd(nanmean(nanmean(evAmp_5ms_1(:,cfg.layerassign{4},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_1(:,cfg.layerassign{4},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_1(:,cfg.layerassign{4},:),2),3));...
    nanstd(nanmean(nanmean(evAmp_5ms_1(:,cfg.layerassign{5},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_1(:,cfg.layerassign{5},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_1(:,cfg.layerassign{5},:),2),3));...
    nanstd(nanmean(nanmean(evAmp_5ms_1(:,cfg.layerassign{6},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_1(:,cfg.layerassign{6},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_1(:,cfg.layerassign{6},:),2),3))];

% Plot
fH = figure();
axes()
set(gca,'ColorOrder',cm([1,3,5],:),'NextPlot','replacechildren')
errorbar(tmp,tmps./(sqrt(9))) % mean +- SEM
set(gca,'YLim',[12 28])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_TFR_evokedAmp_allstimLayer_allDelay_shank1')
close(fH)

%
% Shank 2
%
% Mean
tmp = [nanmean(nanmean(nanmean(evAmp_5ms_2(:,cfg.layerassign{1},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_2(:,cfg.layerassign{1},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_2(:,cfg.layerassign{1},:),2),3));...
    nanmean(nanmean(nanmean(evAmp_5ms_2(:,cfg.layerassign{2},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_2(:,cfg.layerassign{2},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_2(:,cfg.layerassign{2},:),2),3));...
    nanmean(nanmean(nanmean(evAmp_5ms_2(:,cfg.layerassign{3},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_2(:,cfg.layerassign{3},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_2(:,cfg.layerassign{3},:),2),3));...
    nanmean(nanmean(nanmean(evAmp_5ms_2(:,cfg.layerassign{4},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_2(:,cfg.layerassign{4},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_2(:,cfg.layerassign{4},:),2),3));...
    nanmean(nanmean(nanmean(evAmp_5ms_2(:,cfg.layerassign{5},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_2(:,cfg.layerassign{5},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_2(:,cfg.layerassign{5},:),2),3));...
    nanmean(nanmean(nanmean(evAmp_5ms_2(:,cfg.layerassign{6},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_2(:,cfg.layerassign{6},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_2(:,cfg.layerassign{6},:),2),3))];

% Standard deviation
tmps = [nanstd(nanmean(nanmean(evAmp_5ms_2(:,cfg.layerassign{1},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_2(:,cfg.layerassign{1},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_2(:,cfg.layerassign{1},:),2),3));...
    nanstd(nanmean(nanmean(evAmp_5ms_2(:,cfg.layerassign{2},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_2(:,cfg.layerassign{2},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_2(:,cfg.layerassign{2},:),2),3));...
    nanstd(nanmean(nanmean(evAmp_5ms_2(:,cfg.layerassign{3},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_2(:,cfg.layerassign{3},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_2(:,cfg.layerassign{3},:),2),3));...
    nanstd(nanmean(nanmean(evAmp_5ms_2(:,cfg.layerassign{4},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_2(:,cfg.layerassign{4},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_2(:,cfg.layerassign{4},:),2),3));...
    nanstd(nanmean(nanmean(evAmp_5ms_2(:,cfg.layerassign{5},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_2(:,cfg.layerassign{5},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_2(:,cfg.layerassign{5},:),2),3));...
    nanstd(nanmean(nanmean(evAmp_5ms_2(:,cfg.layerassign{6},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_2(:,cfg.layerassign{6},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_2(:,cfg.layerassign{6},:),2),3))];

% Plot
fH = figure();
axes()
set(gca,'ColorOrder',cm([1,3,5],:),'NextPlot','replacechildren')
errorbar(tmp,tmps./(sqrt(9))) % mean +- SEM
set(gca,'YLim',[12 28])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_TFR_evokedAmp_allstimLayer_allDelay_shank2')
close(fH)

%
% ECoG
%
% Mean
tmp = [nanmean(nanmean(nanmean(evAmp_5ms_3(:,cfg.layerassign{1},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_3(:,cfg.layerassign{1},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_3(:,cfg.layerassign{1},:),2),3));...
    nanmean(nanmean(nanmean(evAmp_5ms_3(:,cfg.layerassign{2},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_3(:,cfg.layerassign{2},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_3(:,cfg.layerassign{2},:),2),3));...
    nanmean(nanmean(nanmean(evAmp_5ms_3(:,cfg.layerassign{3},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_3(:,cfg.layerassign{3},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_3(:,cfg.layerassign{3},:),2),3));...
    nanmean(nanmean(nanmean(evAmp_5ms_3(:,cfg.layerassign{4},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_3(:,cfg.layerassign{4},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_3(:,cfg.layerassign{4},:),2),3));...
    nanmean(nanmean(nanmean(evAmp_5ms_3(:,cfg.layerassign{5},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_3(:,cfg.layerassign{5},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_3(:,cfg.layerassign{5},:),2),3));...
    nanmean(nanmean(nanmean(evAmp_5ms_3(:,cfg.layerassign{6},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_3(:,cfg.layerassign{6},:),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_3(:,cfg.layerassign{6},:),2),3))];

% Standard deviation
tmps = [nanstd(nanmean(nanmean(evAmp_5ms_3(:,cfg.layerassign{1},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_3(:,cfg.layerassign{1},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_3(:,cfg.layerassign{1},:),2),3));...
    nanstd(nanmean(nanmean(evAmp_5ms_3(:,cfg.layerassign{2},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_3(:,cfg.layerassign{2},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_3(:,cfg.layerassign{2},:),2),3));...
    nanstd(nanmean(nanmean(evAmp_5ms_3(:,cfg.layerassign{3},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_3(:,cfg.layerassign{3},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_3(:,cfg.layerassign{3},:),2),3));...
    nanstd(nanmean(nanmean(evAmp_5ms_3(:,cfg.layerassign{4},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_3(:,cfg.layerassign{4},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_3(:,cfg.layerassign{4},:),2),3));...
    nanstd(nanmean(nanmean(evAmp_5ms_3(:,cfg.layerassign{5},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_3(:,cfg.layerassign{5},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_3(:,cfg.layerassign{5},:),2),3));...
    nanstd(nanmean(nanmean(evAmp_5ms_3(:,cfg.layerassign{6},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_3(:,cfg.layerassign{6},:),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_3(:,cfg.layerassign{6},:),2),3))];

% Plot
fH = figure();
axes()
set(gca,'ColorOrder',cm([1,3,5],:),'NextPlot','replacechildren')
errorbar(tmp,tmps./(sqrt(7))) % mean +- SEM
set(gca,'YLim',[12 28])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_TFR_evokedAmp_allstimLayer_allDelay_ECoG')
close(fH)

%%
% ---------
% Figure 5E
% ---------
%
% Mean
tmp = [nanmean(nanmean(nanmean(evAmp_5ms_1(:,:,cfg.alpha),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_1(:,:,cfg.alpha),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_1(:,:,cfg.alpha),2),3));...
    nanmean(nanmean(nanmean(evAmp_5ms_1(:,:,cfg.beta),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_1(:,:,cfg.beta),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_1(:,:,cfg.beta),2),3));...
    nanmean(nanmean(nanmean(evAmp_5ms_1(:,:,cfg.logamma),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_1(:,:,cfg.logamma),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_1(:,:,cfg.logamma),2),3));...
    nanmean(nanmean(nanmean(evAmp_5ms_1(:,:,cfg.higamma),2),3)),...
    nanmean(nanmean(nanmean(evAmp_15ms_1(:,:,cfg.higamma),2),3)),...
    nanmean(nanmean(nanmean(evAmp_25ms_1(:,:,cfg.higamma),2),3))];

% Standard deviation
tmps = [nanstd(nanmean(nanmean(evAmp_5ms_1(:,:,cfg.alpha),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_1(:,:,cfg.alpha),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_1(:,:,cfg.alpha),2),3));...
    nanstd(nanmean(nanmean(evAmp_5ms_1(:,:,cfg.beta),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_1(:,:,cfg.beta),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_1(:,:,cfg.beta),2),3));...
    nanstd(nanmean(nanmean(evAmp_5ms_1(:,:,cfg.logamma),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_1(:,:,cfg.logamma),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_1(:,:,cfg.logamma),2),3));...
    nanstd(nanmean(nanmean(evAmp_5ms_1(:,:,cfg.higamma),2),3)),...
    nanstd(nanmean(nanmean(evAmp_15ms_1(:,:,cfg.higamma),2),3)),...
    nanstd(nanmean(nanmean(evAmp_25ms_1(:,:,cfg.higamma),2),3))];

% Plot
fH = figure();
axes()
set(gca,'ColorOrder',cm([1,3,5],:),'NextPlot','replacechildren')
errorbar(tmp,tmps./(sqrt(9))) % mean +- SEM
set(gca,'YLim',[10 30])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_TFR_evokedAmp_allfreqband_allDelay_shank1')
close(fH)

clear tmp tmps

%% Statistics
%
% TWO WAY ANOVA - Evoked Response
%
% Collect data:
% 5 ms
values5 = [(nanmean(nanmean(evAmp_5ms_1(:,cfg.layerassign{1},:),2),3));...
    (nanmean(nanmean(evAmp_5ms_1(:,cfg.layerassign{2},:),2),3));...
    (nanmean(nanmean(evAmp_5ms_1(:,cfg.layerassign{3},:),2),3));...
    (nanmean(nanmean(evAmp_5ms_1(:,cfg.layerassign{4},:),2),3));...
    (nanmean(nanmean(evAmp_5ms_1(:,cfg.layerassign{5},:),2),3));...
    (nanmean(nanmean(evAmp_5ms_1(:,cfg.layerassign{6},:),2),3))];
% 15 ms
values15 = [(nanmean(nanmean(evAmp_15ms_1(:,cfg.layerassign{1},:),2),3));...
    (nanmean(nanmean(evAmp_15ms_1(:,cfg.layerassign{2},:),2),3));...
    (nanmean(nanmean(evAmp_15ms_1(:,cfg.layerassign{3},:),2),3));...
    (nanmean(nanmean(evAmp_15ms_1(:,cfg.layerassign{4},:),2),3));...
    (nanmean(nanmean(evAmp_15ms_1(:,cfg.layerassign{5},:),2),3));...
    (nanmean(nanmean(evAmp_15ms_1(:,cfg.layerassign{6},:),2),3))];
% 25 ms
values25 = [(nanmean(nanmean(evAmp_25ms_1(:,cfg.layerassign{1},:),2),3));...
    (nanmean(nanmean(evAmp_25ms_1(:,cfg.layerassign{2},:),2),3));...
    (nanmean(nanmean(evAmp_25ms_1(:,cfg.layerassign{3},:),2),3));...
    (nanmean(nanmean(evAmp_25ms_1(:,cfg.layerassign{4},:),2),3));...
    (nanmean(nanmean(evAmp_25ms_1(:,cfg.layerassign{5},:),2),3));...
    (nanmean(nanmean(evAmp_25ms_1(:,cfg.layerassign{6},:),2),3))];

% Concatenate all values
values = [values5;values15;values25];
clear values5 values15 values25

% Factor definitions
factor1 = repmat(kron((1:6)',ones(9,1)),3,1);
factor2 = [ones(9*6,1)*5;ones(9*6,1)*15;ones(9*6,1)*25];
grouping = {factor1,factor2};
clear factor1 factor2

% Two-way ANOVA:
[p,tbl,stats] = anovan(values,grouping,'model','interaction','random',2,'varnames',{'Layer','Delay'}); %#ok values used for manuscript text

% Tukey-HSD post-hoc test:
c = multcompare(stats,'Dimension',2); %#ok values used for manuscript text, not further processing

clear p tbl stats c values grouping

%% Split experiments into "super-additive" and "neutral"

% Load data (again)
load('C:\mbv\temp\tfr\res_5ms_shank1','res_5ms')

% Get induced repsonse amplitude
y = squeeze(nanmean(nanmean(max(max(res_5ms(:,:,cfg.supgran,:,cfg.InducedResponse),[],5),[],4),2),3));
clear res_5ms

% k-means clustering with k = 2
[clust_idx,clust_cent] = kmeans(y,2);

% ---------
% Figure 5F
% ---------
%
colmap = [1 0 0; 0 0 1];
fH = figure();
for eID = 1:length(cfg.experiments)
    plot(1,y(eID),'o','MarkerEdgeColor','w','MarkerFaceColor',colmap(clust_idx(eID),:))
    hold on
end
plot([0.5 1.5],[clust_cent(1) clust_cent(1)],'--r')
plot([0.5 1.5],[clust_cent(2) clust_cent(2)],'--b')
set(gca,'YLim',[0 20],'XLim',[0 2])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_split_experiments_res_5ms')
close(fH)

% Define super-additive and neutral experiments according to clustering:
cfg.SupraAddExps = clust_idx == 1;
cfg.NonRespoExps = clust_idx == 2;

%
% Acoustic induced response for comparison
%
% Load data
[pow_click,~] = Load_TFR_Acoustic(cfg);

% Get induced response amplitude
y2 = squeeze(mean(max(max(pow_click(:,cfg.supgran,:,cfg.InducedResponse),[],4),[],3),2));
clear pow_click

% Plot
fH = figure();
for eID = 1:length(cfg.experiments)
    plot(1,y2(eID),'o','MarkerEdgeColor','w','MarkerFaceColor',colmap(clust_idx(eID),:))
    hold on
end
plot([0.5 1.5],[mean(y2(clust_idx==1)) mean(y2(clust_idx==1))],'--r')
plot([0.5 1.5],[mean(y2(clust_idx==2)) mean(y2(clust_idx==2))],'--b')
set(gca,'YLim',[0 20],'XLim',[0 2])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_split_experiments_click_induced')
close(fH)

clear y clust_idx clust_cent eID y2

%% Difference data pictures
% ---------
% Figure 6A
% ---------
%
% Load 5ms difference data
load('C:\mbv\temp\tfr\res_5ms_shank1','res_5ms')
temp = res_5ms; clear res_5ms
load('C:\mbv\temp\tfr\res_5ms_shank2','res_5ms')
temp = cat(3,temp,res_5ms); clear res_5ms
res_5ms = temp;
clear temp

% Plot
for stimL = 1:6
    % Shank 1
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{stimL},cfg.supgran,:,:),3),2),1)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800',['C:\mbv\temp\tfr\plots\Diff_5ms_SA_Shank1_stimLayer' num2str(stimL) '_pow_supragranular'])
    close(fH)
    % Shank 2
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{stimL},cfg.supgran+16,:,:),3),2),1)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800',['C:\mbv\temp\tfr\plots\Diff_5ms_SA_Shank2_stimLayer' num2str(stimL) '_pow_supragranular'])
    close(fH)
end

clear res_5ms

% ---------
% Figure 6B
% ---------
%
% Load 15ms difference data
load('C:\mbv\temp\tfr\res_15ms_shank1','res_15ms')
temp = res_15ms; clear res_15ms %#ok variable loaded from file
load('C:\mbv\temp\tfr\res_15ms_shank2','res_15ms')
temp = cat(3,temp,res_15ms); clear res_15ms
res_15ms = temp;
clear temp

% Plot
for stimL = 1:6
    % Shank 1
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{stimL},cfg.supgran,:,:),3),2),1)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800',['C:\mbv\temp\tfr\plots\Diff_15ms_SA_Shank1_stimLayer' num2str(stimL) '_pow_supragranular2'])
    close(fH)
    % Shank 2
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{stimL},cfg.supgran+16,:,:),3),2),1)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800',['C:\mbv\temp\tfr\plots\Diff_15ms_SA_Shank2_stimLayer' num2str(stimL) '_pow_supragranular2'])
    close(fH)
end

clear res_15ms

% ---------
% Figure 6C
% ---------
%
% Load 25ms difference data
load('C:\mbv\temp\tfr\res_25ms_shank1','res_25ms')
temp = res_25ms; clear res_25ms %#ok variable loaded from file
load('C:\mbv\temp\tfr\res_25ms_shank2','res_25ms')
temp = cat(3,temp,res_25ms); clear res_25ms
res_25ms = temp;
clear temp

% Plot
for stimL = 1:6
    % Shank 1
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{stimL},cfg.supgran,:,:),3),2),1)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800',['C:\mbv\temp\tfr\plots\Diff_25ms_SA_Shank1_stimLayer' num2str(stimL) '_pow_supragranular2'])
    close(fH)
    % Shank 2
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{stimL},cfg.supgran+16,:,:),3),2),1)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800',['C:\mbv\temp\tfr\plots\Diff_25ms_SA_Shank2_stimLayer' num2str(stimL) '_pow_supragranular2'])
    close(fH)
end

clear res_25ms

%%
% Shank 1
load('C:\mbv\temp\tfr\res_5ms_shank1','res_5ms')
load('C:\mbv\temp\tfr\res_15ms_shank1','res_15ms')
load('C:\mbv\temp\tfr\res_25ms_shank1','res_25ms')

% Initialize variables
evAmp_5ms = zeros(length(cfg.experiments),length(cfg.fi),45);
inAmp_5ms = zeros(length(cfg.experiments),length(cfg.fi),45);
evAmp_15ms = zeros(length(cfg.experiments),length(cfg.fi),45);
inAmp_15ms = zeros(length(cfg.experiments),length(cfg.fi),45);
evAmp_25ms = zeros(length(cfg.experiments),length(cfg.fi),45);
inAmp_25ms = zeros(length(cfg.experiments),length(cfg.fi),45);

% Get amplitudes
for eID = 1:length(cfg.experiments)
    for f = 1:16
        evAmp_5ms(eID,f,:) = nanmean(max(res_5ms(eID,f,cfg.supgran,:,cfg.EvokedResponse),[],5),3);
        inAmp_5ms(eID,f,:) = nanmean(max(res_5ms(eID,f,cfg.supgran,:,cfg.InducedResponse),[],5),3);
        evAmp_15ms(eID,f,:) = nanmean(max(res_15ms(eID,f,cfg.supgran,:,cfg.EvokedResponse),[],5),3);
        inAmp_15ms(eID,f,:) = nanmean(max(res_15ms(eID,f,cfg.supgran,:,cfg.InducedResponse),[],5),3);
        evAmp_25ms(eID,f,:) = nanmean(max(res_25ms(eID,f,cfg.supgran,:,cfg.EvokedResponse),[],5),3);
        inAmp_25ms(eID,f,:) = nanmean(max(res_25ms(eID,f,cfg.supgran,:,cfg.InducedResponse),[],5),3);
    end
end

% ---------
% Figure 6D
% ---------
% 
% Induced amplitude as a function of delta_t
fH = figure();
plot(ones(1,4)*1,nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,:),2),3),'or')
hold on
plot([0.5 1.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,:),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,:),2),3))],':r')
plot(ones(1,5)*1,nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,:),2),3),'ok')
plot([0.5 1.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,:),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,:),2),3))],':k')
plot(ones(1,4)*2,nanmean(nanmean(inAmp_15ms(cfg.SupraAddExps,:,:),2),3),'or')
plot([1.5 2.5],[mean(nanmean(nanmean(inAmp_15ms(cfg.SupraAddExps,:,:),2),3)) mean(nanmean(nanmean(inAmp_15ms(cfg.SupraAddExps,:,:),2),3))],':r')
plot(ones(1,5)*2,nanmean(nanmean(inAmp_15ms(cfg.NonRespoExps,:,:),2),3),'ok')
plot([1.5 2.5],[mean(nanmean(nanmean(inAmp_15ms(cfg.NonRespoExps,:,:),2),3)) mean(nanmean(nanmean(inAmp_15ms(cfg.NonRespoExps,:,:),2),3))],':k')
plot(ones(1,4)*3,nanmean(nanmean(inAmp_25ms(cfg.SupraAddExps,:,:),2),3),'or')
plot([2.5 3.5],[mean(nanmean(nanmean(inAmp_25ms(cfg.SupraAddExps,:,:),2),3)) mean(nanmean(nanmean(inAmp_25ms(cfg.SupraAddExps,:,:),2),3))],':r')
plot(ones(1,5)*3,nanmean(nanmean(inAmp_25ms(cfg.NonRespoExps,:,:),2),3),'ok')
plot([2.5 3.5],[mean(nanmean(nanmean(inAmp_25ms(cfg.NonRespoExps,:,:),2),3)) mean(nanmean(nanmean(inAmp_25ms(cfg.NonRespoExps,:,:),2),3))],':k')
plot([0 4],[0 0],'--k')
set(gca,'YLim',[-5 12],'XLim',[0 4])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_diff_ampvs0_freqmean_alldelay_induced')
close(fH)

% ---------
% Figure 6E
% ---------
% 
% Induced amplitude as a function of stimulated layer - 5 ms data
fH = figure();
% Layer 1
plot(ones(1,4)*1,nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{1},:),2),3),'or')
hold on
plot([0.5 1.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{1},:),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{1},:),2),3))],':r')
plot(ones(1,5)*1,nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{1},:),2),3),'ok')
plot([0.5 1.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{1},:),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{1},:),2),3))],':k')
% Layer 2
plot(ones(1,4)*2,nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{2},:),2),3),'or')
plot([1.5 2.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{2},:),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{2},:),2),3))],':r')
plot(ones(1,5)*2,nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{2},:),2),3),'ok')
plot([1.5 2.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{2},:),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{2},:),2),3))],':k')
% Layer 3
plot(ones(1,4)*3,nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{3},:),2),3),'or')
plot([2.5 3.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{3},:),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{3},:),2),3))],':r')
plot(ones(1,5)*3,nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{3},:),2),3),'ok')
plot([2.5 3.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{3},:),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{3},:),2),3))],':k')
% Layer 4
plot(ones(1,4)*4,nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{4},:),2),3),'or')
plot([3.5 4.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{4},:),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{4},:),2),3))],':r')
plot(ones(1,5)*4,nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{4},:),2),3),'ok')
plot([3.5 4.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{4},:),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{4},:),2),3))],':k')
% Layer 5
plot(ones(1,4)*5,nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{5},:),2),3),'or')
plot([4.5 5.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{5},:),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{5},:),2),3))],':r')
plot(ones(1,5)*5,nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{5},:),2),3),'ok')
plot([4.5 5.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{5},:),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{5},:),2),3))],':k')
% Layer 6
plot(ones(1,4)*6,nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{6},:),2),3),'or')
plot([5.5 6.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{6},:),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{6},:),2),3))],':r')
plot(ones(1,5)*6,nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{6},:),2),3),'ok')
plot([5.5 6.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{6},:),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{6},:),2),3))],':k')
set(gca,'YLim',[-5 12],'XLim',[0 7])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_diff_ampvs0_freqmean_alldlayers_induced')
close(fH)

% ---------
% Figure 6F
% ---------
% 
% Induced amplitude as a function of frequency band - 5 ms data
fH = figure();
plot(ones(1,4)*1,nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,cfg.alpha),2),3),'or')
hold on
plot([0.5 1.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,cfg.alpha),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,cfg.alpha),2),3))],':r')
plot(ones(1,4)*2,nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,cfg.beta),2),3),'or')
plot([1.5 2.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,cfg.beta),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,cfg.beta),2),3))],':r')
plot(ones(1,4)*3,nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,cfg.logamma),2),3),'or')
plot([2.5 3.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,cfg.logamma),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,cfg.logamma),2),3))],':r')
plot(ones(1,4)*4,nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,cfg.higamma),2),3),'or')
plot([3.5 4.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,cfg.higamma),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,cfg.higamma),2),3))],':r')
plot(ones(1,5)*1,nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,cfg.alpha),2),3),'ok')
plot([0.5 1.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,cfg.alpha),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,cfg.alpha),2),3))],':k')
plot(ones(1,5)*2,nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,cfg.beta),2),3),'ok')
plot([1.5 2.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,cfg.beta),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,cfg.beta),2),3))],':k')
plot(ones(1,5)*3,nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,cfg.logamma),2),3),'ok')
plot([2.5 3.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,cfg.logamma),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,cfg.logamma),2),3))],':k')
plot(ones(1,5)*4,nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,cfg.higamma),2),3),'ok')
plot([3.5 4.5],[mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,cfg.higamma),2),3)) mean(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,cfg.higamma),2),3))],':k')
plot([0 5],[0 0],'--k')
set(gca,'YLim',[-5 12],'XLim',[0 5])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_diff_ampvsfreq_5ms_induced')
close(fH)

%% Statistics
clear p ci stats

% Perform one-sample t-tests to see if the mean of each group differs 
% significantly from 0

% Factor: delta_t
[~,p(1),ci(1,:),stats(1)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,:),2),3));
[~,p(2),ci(2,:),stats(2)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,:),2),3));
[~,p(3),ci(3,:),stats(3)] = ttest(nanmean(nanmean(inAmp_15ms(cfg.SupraAddExps,:,:),2),3));
[~,p(4),ci(4,:),stats(4)] = ttest(nanmean(nanmean(inAmp_15ms(cfg.NonRespoExps,:,:),2),3));
[~,p(5),ci(5,:),stats(5)] = ttest(nanmean(nanmean(inAmp_25ms(cfg.SupraAddExps,:,:),2),3));
[~,p(6),ci(6,:),stats(6)] = ttest(nanmean(nanmean(inAmp_25ms(cfg.NonRespoExps,:,:),2),3));

% Factor: layer
[~,p(7),ci(1,:),stats(7)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{1},:),2),3));
[~,p(8),ci(2,:),stats(8)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{1},:),2),3));
[~,p(9),ci(3,:),stats(9)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{2},:),2),3));
[~,p(10),ci(4,:),stats(10)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{2},:),2),3));
[~,p(11),ci(5,:),stats(11)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{3},:),2),3));
[~,p(12),ci(6,:),stats(12)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{3},:),2),3));
[~,p(13),ci(7,:),stats(13)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{4},:),2),3));
[~,p(14),ci(8,:),stats(14)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{4},:),2),3));
[~,p(15),ci(9,:),stats(15)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{5},:),2),3));
[~,p(16),ci(10,:),stats(16)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{5},:),2),3));
[~,p(17),ci(11,:),stats(17)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,cfg.layerassign{6},:),2),3));
[~,p(18),ci(12,:),stats(18)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,cfg.layerassign{6},:),2),3));

% Factor: frequency band
[~,p(19),ci(19,:),stats(19)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,cfg.alpha),2),3));
[~,p(20),ci(20,:),stats(20)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,cfg.alpha),2),3));
[~,p(21),ci(21,:),stats(21)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,cfg.beta),2),3));
[~,p(22),ci(22,:),stats(22)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,cfg.beta),2),3));
[~,p(23),ci(23,:),stats(23)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,cfg.logamma),2),3));
[~,p(24),ci(24,:),stats(24)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,cfg.logamma),2),3));
[~,p(25),ci(25,:),stats(25)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.SupraAddExps,:,cfg.higamma),2),3));
[~,p(26),ci(26,:),stats(26)] = ttest(nanmean(nanmean(inAmp_5ms(cfg.NonRespoExps,:,cfg.higamma),2),3)); %#ok yes, ci and stats are unused...

% Benjamini-Hochberg correction for multiple testing
psorted = sort(p);
pcorrected = zeros(length(psorted),1);
for idx = 1:length(psorted)
    pk = idx./26.*0.05;
    pcorrected(idx) = psorted(idx) <= pk;
end

fprintf('%5.4f\n',p(:))
clear psorted pcorrected p pk ci stats

%% 3-way ANOVA of induced power in experiments classified as super-additive
%
% Factors :
%   - Stimulation depth (Layer)
%   - Delay between A & E stimulation
%   - Frequency band

factor1 = repmat((kron((1:6)',ones(4*4,1))),3,1);
factor2 = [ones(4*4*6,1)*5;ones(4*4*6,1)*15;ones(4*4*6,1)*25];
factor3 = repmat({'alpha','beta','logamma','higamma'}',4*6*3,1);
grouping = {factor1,factor2,factor3};
clear factor1 factor2 factor3

values_5ms = [max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{1},cfg.supgran,cfg.alpha,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{1},cfg.supgran,cfg.beta,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{1},cfg.supgran,cfg.logamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{1},cfg.supgran,cfg.higamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{2},cfg.supgran,cfg.alpha,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{2},cfg.supgran,cfg.beta,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{2},cfg.supgran,cfg.logamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{2},cfg.supgran,cfg.higamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{3},cfg.supgran,cfg.alpha,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{3},cfg.supgran,cfg.beta,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{3},cfg.supgran,cfg.logamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{3},cfg.supgran,cfg.higamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{4},cfg.supgran,cfg.alpha,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{4},cfg.supgran,cfg.beta,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{4},cfg.supgran,cfg.logamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{4},cfg.supgran,cfg.higamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{5},cfg.supgran,cfg.alpha,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{5},cfg.supgran,cfg.beta,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{5},cfg.supgran,cfg.logamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{5},cfg.supgran,cfg.higamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{6},cfg.supgran,cfg.alpha,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{6},cfg.supgran,cfg.beta,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{6},cfg.supgran,cfg.logamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_5ms(cfg.SupraAddExps,cfg.layerassign{6},cfg.supgran,cfg.higamma,600:800),3),2),[],4),[],5)];

values_15ms = [max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{1},cfg.supgran,cfg.alpha,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{1},cfg.supgran,cfg.beta,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{1},cfg.supgran,cfg.logamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{1},cfg.supgran,cfg.higamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{2},cfg.supgran,cfg.alpha,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{2},cfg.supgran,cfg.beta,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{2},cfg.supgran,cfg.logamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{2},cfg.supgran,cfg.higamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{3},cfg.supgran,cfg.alpha,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{3},cfg.supgran,cfg.beta,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{3},cfg.supgran,cfg.logamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{3},cfg.supgran,cfg.higamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{4},cfg.supgran,cfg.alpha,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{4},cfg.supgran,cfg.beta,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{4},cfg.supgran,cfg.logamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{4},cfg.supgran,cfg.higamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{5},cfg.supgran,cfg.alpha,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{5},cfg.supgran,cfg.beta,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{5},cfg.supgran,cfg.logamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{5},cfg.supgran,cfg.higamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{6},cfg.supgran,cfg.alpha,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{6},cfg.supgran,cfg.beta,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{6},cfg.supgran,cfg.logamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_15ms(cfg.SupraAddExps,cfg.layerassign{6},cfg.supgran,cfg.higamma,600:800),3),2),[],4),[],5)];

values_25ms = [max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{1},cfg.supgran,cfg.alpha,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{1},cfg.supgran,cfg.beta,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{1},cfg.supgran,cfg.logamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{1},cfg.supgran,cfg.higamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{2},cfg.supgran,cfg.alpha,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{2},cfg.supgran,cfg.beta,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{2},cfg.supgran,cfg.logamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{2},cfg.supgran,cfg.higamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{3},cfg.supgran,cfg.alpha,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{3},cfg.supgran,cfg.beta,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{3},cfg.supgran,cfg.logamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{3},cfg.supgran,cfg.higamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{4},cfg.supgran,cfg.alpha,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{4},cfg.supgran,cfg.beta,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{4},cfg.supgran,cfg.logamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{4},cfg.supgran,cfg.higamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{5},cfg.supgran,cfg.alpha,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{5},cfg.supgran,cfg.beta,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{5},cfg.supgran,cfg.logamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{5},cfg.supgran,cfg.higamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{6},cfg.supgran,cfg.alpha,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{6},cfg.supgran,cfg.beta,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{6},cfg.supgran,cfg.logamma,600:800),3),2),[],4),[],5);...
    max(max(nanmean(nanmean(res_25ms(cfg.SupraAddExps,cfg.layerassign{6},cfg.supgran,cfg.higamma,600:800),3),2),[],4),[],5)];

values = [values_5ms;values_15ms;values_25ms];

[p,tbl,stats] = anovan(values,grouping,'model','full','random',2,'varnames',{'Layer','Delay','FreqBand'}); %#ok the values here are used for the manuscript text but not in any further calculation
clear values values_5ms values_15ms values_25ms p tbl stats grouping

