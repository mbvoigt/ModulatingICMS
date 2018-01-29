%% Plot MUA - Acoustic
%
% Subfunction for ModulatingICMS_main
%
% Plots multi unit activity responses for acoustically only stimulated
% trials
%
% Plots figure 2A
%
function Plot_MUA_Acoustic(~,MUA_click)

% ---------
% FIGURE 2A - Multi unit activity example
% ---------
for example_exp = 1:length(cfg.experiments)
    
    % Threshold MUA activity to detect spikes
    % (on the first shank only)
    spikematrixlog = spikes(squeeze(MUA_click(example_exp,1:16,:,:)));
    
    % MUA spike raster
    fH = figure();
    hold on
    for chan = 1:16
        for sw = 1:30
            plotdatax = find(spikematrixlog(chan,sw,:));
            plotdatay = -1*(sw+30*chan)*ones(length(find(spikematrixlog(chan,sw,:))),1);
            scatter(gca,plotdatax,plotdatay,'.','MarkerEdgeColor','k','MarkerFaceColor','k')
        end
    end
    for chan = 2:16
        plot([0 22001],[-30*chan -30*chan],':k')
    end
    plot([0 0],[-30 0],'k')
    plot([0 2200],[0 0],'k')
    axis off
    print(fH,'-dpng','-r1800',['C:\mbv\temp\tfr\plots\TEI_click_MUA_raster_' cfg.experiments{example_exp}])
    close(fH)
    
    % MUA PSTH
    % Calculate spike times & create histogram data
    spTS = [];
    for ch=1:size(spikematrixlog,1)
        for sw=1:size(spikematrixlog,2)
            spTS = [spTS; find(spikematrixlog(ch,sw,:))]; %#ok no. of spikes isn't knwown beforehand
        end
    end
    spTS = spTS/22000*1000; % convert samples to ms
    [n,xout] = hist(spTS,size(spikematrixlog,3)/SR*1000); % get histogram data
    
    % Plot
    fH = figure();
    bar(xout,n,'k')
    set(gca,'XTick',0:100:1100,'XTickLabel',-500:100:600)
    print(fH,'-dpdf','-r1800',['C:\mbv\temp\tfr\plots\TEI_click_MUA_PSTH_' cfg.experiments{example_exp}])
    close(fH)
end
% Clear temporary/helper variables
clear example_exp chan sw plotdatax plotdatay spikematrixlog n xout
