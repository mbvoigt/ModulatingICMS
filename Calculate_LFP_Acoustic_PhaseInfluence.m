%% Calculate LFP - Acoustic - Pre-stimulus phase influence
%
% Subfunction for ModulatingICMS_main
%
% Calculates pre-stimulus phase influence of acoustically only 
% stimulated trials
%
% Plots figures 7A, 7B, 7D, 7F
%
function Calculate_LFP_Acoustic_PhaseInfluence(cfg,LFP_click)

% Remove all data after 50 ms post-stimulation to save memory
LFP_click = single(LFP_click(:,:,:,1:cfg.sweepBefore+2200));

% ---------
% FIGURE 7A - Example: LFP colored according to phase
% ---------
% take example data
tmp = squeeze(nanmean(nanmean(LFP_click(:,4,:,11000-1100:end),3),1));
% calculate phase
phtmp = hilbert(tmp);
phtmp = angle(phtmp);
% 8 phase bins from -pi to pi
binrange = -pi():2*pi()/8:pi();
% bin phase values into these 8 bins
[~,~,bin] = histcounts(phtmp,binrange);
% define a colormap with 8 steps
cm = jet(8);

% Plotting
fH = figure();
hold on
for tidx = 1:length(tmp)
plot(tidx,tmp(tidx),'.','Color',cm(bin(tidx),:),'MarkerSize',9)
end
plot([1100 1100],[-250 100],'--k') % 0 ms
plot([1100+110 1100+110],[-250 100],'--k') % 5 ms
plot([1100+330 1100+330],[-250 100],'--k') % 15 ms
plot([1100+550 1100+550],[-250 100],'--k') % 25 ms
plot([0 1100],[-250 -250],'k','LineWidth',2) % 50 ms
plot([0 0],[-250 -150],'k','LineWidth',2) % 100 µV
colormap(cm)
colorbar()
fH.Renderer = 'painters';
print(fH,'-dpdf','-r300','C:\mbv\temp\tfr\plots\LFP_click_phaseExample_trace')
close(fH)
% ---------

% Calculate LFP phase
phase = zeros(length(cfg.experiments),length(cfg.ch2analyse),30,cfg.sweepBefore+2200);
for eID = 1:length(cfg.experiments)
    for ch = 1:cfg.noChannels
        fprintf('%s | %2d\n',cfg.experiments{eID},ch)
        for sw = 1:30
            phase(eID,ch,sw,:) = hilbert(LFP_click(eID,ch,sw,:));
        end
    end
end
phase = angle(phase);

% Z-scoring of the LFP amplitude
zscore = zeros(length(cfg.experiments),length(cfg.ch2analyse),30,size(LFP_click,4));
for eID = 1:length(cfg.experiments)
    for ch = 1:cfg.noChannels
        fprintf('%s | %2d\n',cfg.experiments{eID},ch)
        for sw = 1:30
            zscore(eID,ch,sw,:) = (squeeze(LFP_click(eID,ch,sw,:))-(squeeze(mean(LFP_click(eID,ch,sw,1:cfg.sweepBefore),4))))./(squeeze(std(LFP_click(eID,ch,sw,1:cfg.sweepBefore),[],4)));
        end
    end
end

% Determine single-trial response amplitude
amp = zeros(length(cfg.experiments),length(cfg.ch2analyse),30);
for eID = 1:length(cfg.experiments)
    for ch = 1:cfg.noChannels
        fprintf('%s | %2d\n',cfg.experiments{eID},ch)
        for sw = 1:30
            amp(eID,ch,sw) = max(abs(zscore(eID,ch,sw,cfg.sweepBefore:cfg.sweepBefore+50*22)),[],4);
        end
    end
end

%%
% ---------
% FIGURE 7B - Example: Phase influence at delta_phi = 0 ms
% ---------
% Example channel:
ch = 7;

% Time point
t = 500*22;

% Plotting
fH = figure('Name','Example phase influence');
tmpph = squeeze(phase(:,ch,:,t));
tmpamp = squeeze(amp(:,ch,:));
scatter(tmpph(:),tmpamp(:),'ok')
hold on
plot([-pi() pi()],[0 0],'--k')  % 0-line
% average with respect to phase into 8 bins
[~,~,bin] = histcounts(tmpph(:),binrange);
averaged = accumarray(bin, tmpamp(:),[],@mean); %  average per bin
averagedstd = accumarray(bin, tmpamp(:),[],@std);   % standard deviation per bin
errorbar(0.5*(binrange(1:end-1)+binrange(2:end)),averaged,averagedstd,'k','LineWidth',2)
set(gca,'XTick',-pi():2*pi()/8:pi())
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\LFP_click_phaseScatter_0ms')
close(fH)

% ---------
% FIGURE 7B - Example: Phase influence at delta_phi = -300 ms
% ---------
% Time point
t = 200*22;

% Plotting
fH = figure('Name','Example phase influence -300ms');
tmpph = squeeze(phase(:,ch,:,t));
scatter(tmpph(:),tmpamp(:),'ok')
hold on
plot([-pi() pi()],[0 0],'--k')
% average with respect to phase into 8 bins
[~,~,bin] = histcounts(tmpph(:),binrange);
tmpamp(bin==0)=[];
bin(bin==0) = [];
averaged = accumarray(bin, tmpamp(:),[],@mean); %  average per bin
averagedstd = accumarray(bin, tmpamp(:),[],@std);   % standard deviation per bin
errorbar(0.5*(binrange(1:end-1)+binrange(2:end)),averaged,averagedstd,'k','LineWidth',2)
set(gca,'XTick',-pi():2*pi()/8:pi())
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\LFP_click_phaseScatter_-300ms')
close(fH)
% ---------

% Binning - Statistics
% Calculating the Modulation Index (MI)
% as range (max - min) of the phase bin means
tpoints=100:25:500;
binrange = -pi():2*pi()/8:pi();
MI = NaN(length(tpoints),cfg.noChannels);
for tidx = 1:length(tpoints)
    t = tpoints(tidx)*22;
    for ch = 1:cfg.noChannels
        tmpph = squeeze(phase(:,ch,:,t));
        tmpamp = squeeze(amp(:,ch,:));
        [~,~,bin] = histcounts(tmpph(:),binrange);
        tmpamp(bin==0)=[];
        bin(bin==0) = [];
        averaged = accumarray(bin, tmpamp(:),[],@mean);
        if ~isempty(averaged)
            MI(tidx,ch) = range(averaged);
        end
    end
end

% ---------
% FIGURE 7D - Modulation Index against delta_phi
% ---------
fH = figure();
boxplot(MI(:,1:16)','Color','k')    % MI of first 16 electrodes (= shank 1)
set(gca,'YLim',[0 9])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\LFP_click_ModulationIndex_Shank1')
close(fH)
% ---------

% ---------
% FIGURE 7F - Modulation Index against recorded cortical layer
% ---------
% t = -300 ms
fH = figure();
plot([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6],MI(5,1:16),'o','MarkerFaceColor',cfg.col_acoustic,'MarkerEdgeColor','none') % Shank 1
hold on
plot([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6],MI(5,17:32),'o','MarkerFaceColor',cfg.col_acousticdark,'MarkerEdgeColor','none') % Shank 2
set(gca,'XLim',[0 7],'YLim',[0 9])
plot(accumarray([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6]', squeeze(MI(5,1:16)),[],@mean),'Color',cfg.col_acoustic);
plot(accumarray([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6]', squeeze(MI(5,17:32)),[],@mean),'Color',cfg.col_acousticdark);
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\LFP_click_ModulationIndex_Layers_-300ms')
close(fH)

% t = 0 ms
fH = figure();
plot([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6],MI(17,1:16),'o','MarkerFaceColor',cfg.col_acoustic,'MarkerEdgeColor','none') % Shank 1
hold on
plot([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6],MI(17,17:32),'o','MarkerFaceColor',cfg.col_acousticdark,'MarkerEdgeColor','none') % Shank 2
set(gca,'XLim',[0 7],'YLim',[0 9])
plot(accumarray([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6]', squeeze(MI(17,1:16)),[],@mean),'Color',cfg.col_acoustic);
plot(accumarray([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6]', squeeze(MI(17,17:32)),[],@mean),'Color',cfg.col_acousticdark);
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\LFP_click_ModulationIndex_Layers_0ms')
close(fH)
% ---------
