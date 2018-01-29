%% Plot LFP - AEES
%
% Subfunction for ModulatingICMS_main
%
% Plots local field potential responses for electrically stimulated
% trials and analyzes surface potential in response to varying ICMS
% stimulation depths
%
% Plots figures 4G, 4H, 4I, 4J, 4K 
%
function Plot_LFP_AEES(cfg,LFP_AEES)

% Experiments with ID 1 & 2 had channel 33 as the current monitor and no
% ECoG grid recording
% -> Remove all ECoG channels in these experiments!
LFP_AEES([1,2],:,33:48,:) = NaN;

% ---------
% FIGURE 4G - Example electric ECoG trace
% ---------
fH = figure();
plot(squeeze(nanmean(nanmean(LFP_AEES(:,cfg.layerassign{2},3+32,450*22:700*22),2),1)),'k')
hold on
plot([0 50*22],[-100 -100],'k','LineWidth',2)
plot([0 0],[-100 -50],'k','LineWidth',2)
plot([50*22 50*22],[-100 100],':k')
plot([90*22 90*22],[-100 100],':k')
plot([160*22 160*22],[-100 100],':k')
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_AEES_ECoG_exampletrace')
close(fH)


% ---------
% FIGURE 4H - Electric ECoG data
% ---------
% Plot surface data on a grid
[X1,X2] = ndgrid(1:4);
[Y1,Y2] = ndgrid(1:0.5:4);
% Channel mapping to Blackrock grid electrode
idx = [13,14,15,16,12,11,10,9,4,3,2,1,5,6,7,8]+32;

% Fine time resolution for AEES
t = 495*22:5*22:540*22;

% Plot for ICMS in all 6 layers
for f = 1:6
    for id = 1:length(t)
        fH=figure();
        Z = reshape(nanmean(nanmean(LFP_AEES(:,cfg.layerassign{f},idx,t(id)),1),2),4,4);
        F = griddedInterpolant(X1,X2,Z);
        V = F(Y1,Y2);
        pcolor(Y1,Y2,V)
        shading interp
        colormap(jet)
        caxis([-100 100])
        axis off
        print(fH,'-dpng','-r1800',['C:\mbv\temp\tfr\plots\TEI_AEES_ECoG_fine_stimLayer' num2str(f) '_' num2str(t(id)/22)])
        close(fH)
    end
end
% Example chosen for Figure 4G: Stimulation in Layer II

% ---------
% FIGURE 4I - ECoG peak amplitude for stimulation in different layers
% ---------
for f = 1:6
    fH=figure();
    Z = reshape(max(nanmean(nanmean(LFP_AEES(:,cfg.layerassign{f},idx,500*22:550*22),1),2),[],4),4,4);
    F = griddedInterpolant(X1,X2,Z);
    V = F(Y1,Y2);
    pcolor(Y1,Y2,V)
    shading interp
    colormap(jet)
    caxis([-100 100])
    axis off
    print(fH,'-dpng','-r1800',['C:\mbv\temp\tfr\plots\TEI_AEES_ECoG_grandMean_stimLayer' num2str(f) '_max'])
    close(fH)
end

% ---------
% FIGURE 4J - Bar chart of ECoG peak amplitudes
% ---------
% Initialize variable
amp = zeros(length(cfg.experiments),6,16);
% Calculate evoked amplitude
for f = 1:6
    amp(:,f,:) = squeeze(max(nanmean(LFP_AEES(:,cfg.layerassign{f},33:48,500*22:550*22),2),[],4));
end

% Take the mean of all 16 ECoG electrodes
tmp = [nanmean(squeeze(amp(:,1,:)),2),...
    nanmean(squeeze(amp(:,2,:)),2),...
    nanmean(squeeze(amp(:,3,:)),2),...
    nanmean(squeeze(amp(:,4,:)),2),...
    nanmean(squeeze(amp(:,5,:)),2),...
    nanmean(squeeze(amp(:,6,:)),2)];
% Plot
fH = figure();
bar(nanmean(tmp),'FaceColor',cfg.col_electric,'EdgeColor',cfg.col_electricdark)
hold on
errorbar(nanmean(tmp),nanstd(tmp)/(sqrt(7)),'Color',cfg.col_electricdark) % Mean of experiments +- SEM
set(gca,'YLim',[0 60])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_AEES_ECoG_LFP_evokedampvslayer')
close(fH)

% Statistics: 1-way ANOVA
[p,tbl,stats] = anova1(tmp); %#ok Results used for manuscript text but not for further calculations


% ---------
% FIGURE 4K - Model for ECoG amplitude in response to ICMS of varying depth
% ---------
% x = depth below pia in µm
x = 0:2500;
% 1st component: linear drop of weight
lindrop = linspace(1,0,length(x));
% 2nd component: normalized gaussian distribution with µ = 1000 and sigma = 400
m = 1000;
s = 500;
gaus = (1/(sqrt(2.*pi().*s^2))).*exp(-((x-m).^2/(2.*s^2)));
gaus = gaus./(1*max(gaus));

% Plot
fH = figure();
plot(x,gaus,'--k')
hold on
plot(x,lindrop,':k')
plot(x,gaus+lindrop,'k')
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_AEES_ECoG_amplitudemodel')
close(fH)







