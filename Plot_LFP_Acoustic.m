%% Plot LFP - Acoustic
%
% Subfunction for ModulatingICMS_main
%
% Plots local field potential responses for acoustically stimulated
% trials
%
% Plots figures 2B, 3H, 3I 
%
function Plot_LFP_Acoustic(~,LFP_click)

% ---------
% FIGURE 2B - Local field potential example
% ---------
for example_exp = 1:length(cfg.experiments)
    
    % Local field potential data 
    % of first shank
    LFP = squeeze(mean(LFP_click(example_exp,1:16,:,:),3));
    
    fH = figure();
    fH.Renderer = 'Painters';
    hold on
    for chan = 1:16
        plot(squeeze(LFP(chan,:))-(200*chan),'k')
    end
    plot([100 100],[-500 0],'k','LineWidth',2)
    set(gca,'XLim',[0 22001],'XTick',0:2200:22001,'XTickLabel',-500:100:550,'YLim',[(-17*200) 0])
    print(fH,'-dpdf','-r1800',['C:\mbv\temp\tfr\plots\TEI_click_LFP_' cfg.experiments{example_exp}])
    close(fH)
    
end
% Clear temporary/helper variables
clear example_exp chan fH

% ---------
% FIGURE 3H - Example acoustic ECoG trace
% ---------
fH = figure();
plot(squeeze(nanmean(nanmean(LFP_click(:,2+32,:,450*22:700*22),3),1)),'k')
hold on
plot([0 50*22],[-150 -150],'k','LineWidth',2)
plot([0 0],[-150 -50],'k','LineWidth',2)
plot([50*22 50*22],[-150 150],':k')
plot([90*22 90*22],[-150 150],':k')
plot([160*22 160*22],[-150 150],':k')
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\TEI_click_ECoG_exampletrace')
close(fH)

% ---------
% FIGURE 3I - Acoustic ECoG data
% ---------
% Plot surface data on a grid
[X1,X2] = ndgrid(1:4);
[Y1,Y2] = ndgrid(1:0.5:4);
% Channel mapping to Blackrock grid electrode
idx = [13,14,15,16,12,11,10,9,4,3,2,1,5,6,7,8]+32;

%{
%%%
% Code to check the correct orientation of plots
%%%
LFP_t = zeros(1,48,1,10);
LFP_t(1,1+32,1,:) = 1;
LFP_t(1,6+32,1,:) = 1;
pcolor(X1,X2,reshape(nanmean(nanmean(LFP_t(:,idx,:,10),3),1),4,4))
shading interp
%%%
%}

% Coarse time resolution
t = 500*22:10*22:740*22;

for id = 1:length(t)
    fH = figure();
    Z = reshape(nanmean(nanmean(LFP_click(:,idx,:,t(id)),3),1),4,4);
    F = griddedInterpolant(X1,X2,Z);
    V = F(Y1,Y2);
    pcolor(Y1,Y2,V);
    shading interp
    colormap(jet)
    caxis([-100 100])
    axis off
    print(fH,'-dpng','-r1800',['C:\mbv\temp\tfr\plots\TEI_click_ECoG_coarse_' num2str(t(id)/22)])
    close(fH)
end

% Fine time resolution
t = 495*22:2*22:544*22;

for id = 1:length(t)
    fH = figure();
    Z = reshape(nanmean(nanmean(LFP_click(:,idx,:,t(id)),3),1),4,4);
    F = griddedInterpolant(X1,X2,Z);
    V = F(Y1,Y2);
    pcolor(Y1,Y2,V);
    shading interp
    colormap(jet)
    caxis([-100 100])
    axis off
    print(fH,'-dpng','-r1800',['C:\mbv\temp\tfr\plots\TEI_click_ECoG_fine_' num2str(t(id)/22)])
    close(fH)
end