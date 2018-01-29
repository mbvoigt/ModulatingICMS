%% Analyse TFR - AEES
%
% Subfunction for ModulatingICMS_main
%
% Analyzed TFRs for electrically stimulated trials & plots further figures
%
% Plots figures 4D, 4E, 4F
%
function Analyze_TFR_AEES(cfg,pow_AEES)

% Colormap for lines showing different layers stimulated
cm = flipud([20/255 11/255 0/255;...
    100/255 55/255 0/255;...
    160/255 88/255 0/255;...
    255/255 140/255 0/255;...
    255/255 160/255 45/255;...
    255/255 205/255 145/255]);

%% Evoked response amplitude
% Initialize variable
evokedAmplitude = NaN(length(cfg.experiments),length(cfg.fi),cfg.noChannels,length(cfg.FTcfg.foi));
% Determine amplitude
for eID = 1:length(cfg.experiments)
    for f = 1:length(cfg.fi)
        for ch = 1: cfg.noChannels
            evokedAmplitude(eID,f,ch,:) = nanmean(pow_AEES(eID,f,ch,:,500:600),5);
        end
    end
end

% ---------
% FIGURE 4E - Evoked response as a function of cortical layer
% ---------
% Shank 1
fH = figure();
for f = 1:6
tmp = [max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},cfg.layerassign{1},:),2),[],4),[],3),...
    max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},cfg.layerassign{2},:),2),[],4),[],3),...
    max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},cfg.layerassign{3},:),2),[],4),[],3),...
    max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},cfg.layerassign{4},:),2),[],4),[],3),...
    max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},cfg.layerassign{5},:),2),[],4),[],3),...
    max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},cfg.layerassign{6},:),2),[],4),[],3)];
errorbar(mean(tmp),std(tmp)/sqrt(9),'Color',cm(f,:)) % mean +- SEM
hold on
end
set(gca,'YLim',[0 13])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_AEES_TFR_evokedAmp_layers_shank1_allstimL')
close(fH)
% ---------

% Shank 2
fH = figure();
for f = 1:6
tmp = [max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},cfg.layerassign{1}+16,:),2),[],4),[],3),...
    max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},cfg.layerassign{2}+16,:),2),[],4),[],3),...
    max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},cfg.layerassign{3}+16,:),2),[],4),[],3),...
    max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},cfg.layerassign{4}+16,:),2),[],4),[],3),...
    max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},cfg.layerassign{5}+16,:),2),[],4),[],3),...
    max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},cfg.layerassign{6}+16,:),2),[],4),[],3)];
errorbar(mean(tmp),std(tmp)/sqrt(9),'Color',cm(f,:)) % mean +- SEM
hold on
end
set(gca,'YLim',[0 13])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_AEES_TFR_evokedAmp_layers_shank2_allstimL')
close(fH)
% ---------

% ---------
% FIGURE 4F - Evoked response as a function of frequency band
% ---------
% Shank 1
fH = figure();
for f = 1:6
tmp = [max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},1:16,cfg.alpha),2),[],4),[],3),...
    max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},1:16,cfg.beta),2),[],4),[],3),...
    max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},1:16,cfg.logamma),2),[],4),[],3),...
    max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},1:16,cfg.higamma),2),[],4),[],3)];
errorbar(mean(tmp),std(tmp)/sqrt(9),'Color',cm(f,:))
hold on
end
set(gca,'YLim',[0 13])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_AEES_TFR_evokedAmp_freqbands_shank1_allstimL')
close(fH)
% ---------

% Shank 2
fH = figure();
for f = 1:6
tmp = [max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},17:32,cfg.alpha),2),[],4),[],3),...
    max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},17:32,cfg.beta),2),[],4),[],3),...
    max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},17:32,cfg.logamma),2),[],4),[],3),...
    max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},17:32,cfg.higamma),2),[],4),[],3)];
errorbar(mean(tmp),std(tmp)/sqrt(9),'Color',cm(f,:))
hold on
end
set(gca,'YLim',[0 13])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_AEES_TFR_evokedAmp_freqbands_shank2_allstimL')
close(fH)
% ---------

% ECoG
fH = figure();
for f = 1:6
tmp = [max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},33:48,cfg.alpha),2),[],4),[],3),...
    max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},33:48,cfg.beta),2),[],4),[],3),...
    max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},33:48,cfg.logamma),2),[],4),[],3),...
    max(max(nanmean(evokedAmplitude(:,cfg.layerassign{f},33:48,cfg.higamma),2),[],4),[],3)];
errorbar(nanmean(tmp),nanstd(tmp)/sqrt(9),'Color',cm(f,:))
hold on
end
set(gca,'YLim',[0 13])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_AEES_TFR_evokedAmp_freqbands_ECoG_allstimL')
close(fH)
% ---------

%%
% ---------
% FIGURE 4D - Single TFR example for 45 µA stimulation
% ---------

% Load raw data
load('E:\ICMS-GP\Data\GP18E15CL\AlphaOmega\MATFiles\ctx_pos1_ss_inf0001.mat','-regexp','^CStim|^CSPK|^CTTL');

% Calculate trigger timestamps
tmstmps1 = (CStimMarker_001(1,:)/(CStimMarker_001_KHz*1000))-CSPK_001_TimeBegin; %#ok Variable is loaded from file
tm = round(tmstmps1*(CSPK_001_KHz*1000));
clear tmstmps1

% Cut out epoochs (=sweeps) of raw data
data_raw = zeros(16,length(tm),1001);
for ch = 1:16
    v = sprintf('CSPK_0%02d',ch);
        for sw = 1:length(tm)
            try
                tmp = eval([v '(tm((1*' num2str(length(tm)) ')-' num2str(length(tm)) '+sw)-cfg.sweepBefore:tm((1*' num2str(length(tm)) ')-' num2str(length(tm)) '+sw)+cfg.sweepAfter)']);
                data_raw(ch,sw,:) = resample(double(tmp),1,22); % Resample to 1kHz!
                clear tmp
            catch
                fprintf(['Channel ' sprintf('%02d',ch) ' | ' sprintf('%02d',sw) ': Couldn''t evaluate sweep!\n'])
            end
        end
end
clear('-regexp','CSPK','CStim','CTTL')
clear ch sw tm v

% Just use sweeps 450:479, corresponding to a current of 45 µA 
data = squeeze(data_raw(:,450:479,:));
clear data_raw

% Blank stimulus artefact (ATTN: range is hardcoded here!)
data = blankstimulus(data,1000,500,503);

% Demean! Subtract mean of baseline from time domain data!
for ch = 1:16
    for sw = 1:30
        baseline = mean(squeeze(data(ch,sw,1:351)));
        data(ch,sw,:) = data(ch,sw,:)-baseline;
    end
end

% Reformat data to fit fieldtrip data format:
% Here data gets reformatted to a structure format,
% at the same time channel 1:16 are now ordered according to the
% NNX electrode order (channelidx)
channelidx = [9,8,10,7,13,4,12,5,15,2,16,1,14,3,11,6];
for tr = 1:30
    for ch=1:16
        dat.trial{tr}(ch,:) = double(data(channelidx(ch),tr,:));
    end
    dat.time{tr} = -0.5:1/1000:0.5;
end
for ch = 1:16
    dat.label{ch}=num2str(ch);
end
dat.fsample=1000;
clear tr ch data channelidx

% Time-Frequency Analysis
TFR = ft_freqanalysis(cfg.FTcfg,dat);
fprintf(repmat('\b',1,826)) % Remove output of freqanalysis
fprintf('\n')
clear dat

% extracting the complex fourier spectrum
fspctrm = single(TFR.fourierspctrm);
clear TFR

% Uncomment for phase-locking factor calculation:
% plfspctrm = fspctrm./abs(fspctrm);
% plfspctrm = squeeze(nanmean(plfspctrm));
% plfspctrm = abs(plfspctrm);
% clear fspctrm

% power spectrum
powspctrm = abs(fspctrm).^2;

% median over trials
powspctrm = squeeze(median(powspctrm,1));

% baseline correction
pow=powspctrm;
clear powspctrm
for ch = 1:16
    for freq = 1:45
        baseline = squeeze(nanmean(pow(ch,freq,cfg.BaselineWindow)));
        pow(ch,freq,:) = 10*log10(pow(ch,freq,:)/baseline);
    end
end
clear baseline ch freq sw

% Plot time-frequency representation examples
% Channel used for Figure 4D: ch = 3

for ch = 1:16
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(pow(ch,:,:)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800',['C:\mbv\temp\tfr\plots\TEI_AEES_highcurrent_pow_' num2str(ch)])
    close(fH)
end

% Uncomment to check phase-locking factor for evoked/induced
% differentiation:
%
% ch = 3;
%
% fH = figure();
% pcolor(-0.5:0.001:0.5,7:2:95,squeeze(plfspctrm(ch,:,:)))
% shading interp
% caxis([0 1])
% colormap(jet)
            