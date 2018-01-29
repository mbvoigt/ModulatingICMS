%% Plot TFR - Acoustic
%
% Subfunction for ModulatingICMS_main
%
% Plots TFRs for acoustically stimulated trials
%
% Plots figures 2C, 3A
%
function Plot_TFR_Acoustic(cfg,data,varargin)
% Plots TFRs
%
% Obligatory input arguments:
%
% data [n-dim array] = data to plot
%
% Name-Value-Pairs:
%
% datatype [string] = 'pow' or 'plf'
% singleexps [logical] = 0 for mean of experiments

%Check inputarguments
options = struct(...
    'datatype','pow');

optionNames = fieldnames(options);
nOptionalArgs = length(varargin);
if round(nOptionalArgs/2)~=nOptionalArgs/2
    error('Please input Name/Value pairs!');
end

for pair = reshape(varargin,2,[])
    inpName = pair{1};
    if any(strcmp(inpName,optionNames))
        options.(inpName) = pair{2};
    else
        error('%s is not a valid parameter name',inpName)
    end
end

% If no experiments are explicitly stated, analyse all experiments
if ~isfield(cfg,'exps2analyse')
    cfg.exps2analyse = 1:size(data,1);
end

if strcmp(options.datatype,'pow')
    % Plotting power spectrum:
    
    % ---------
    % FIGURE 2C - Time-frequency response
    % ---------
    
    % Shank 1
    
    % Layer 1
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(data(cfg.exps2analyse,1:2,:,:),1),2)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\TEI_click_pow_Layer_11')
    close(fH)
    
    % Layer 2
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(data(cfg.exps2analyse,3:4,:,:),1),2)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\TEI_click_pow_Layer_12')
    close(fH)
    
    % Layer 3
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(data(cfg.exps2analyse,5:6,:,:),1),2)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\TEI_click_pow_Layer_13')
    close(fH)
    
    % Layer 4
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(data(cfg.exps2analyse,7:8,:,:),1),2)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\TEI_click_pow_Layer_14')
    close(fH)
    
    % Layer 5
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(data(cfg.exps2analyse,9:12,:,:),1),2)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\TEI_click_pow_Layer_15')
    close(fH)
    
    % Layer 6
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(data(cfg.exps2analyse,13:16,:,:),1),2)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\TEI_click_pow_Layer_16')
    close(fH)
    
    % Grand mean over layers
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(data(cfg.exps2analyse,1:16,:,:),1),2)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\TEI_click_grandMean')
    close(fH)
    
    % Shank 2
    
    % Layer 1
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(data(cfg.exps2analyse,17:18,:,:),1),2)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\TEI_click_pow_Layer_21')
    close(fH)
    
    % Layer 2
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(data(cfg.exps2analyse,19:20,:,:),1),2)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\TEI_click_pow_Layer_22')
    close(fH)
    
    % Layer 3
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(data(cfg.exps2analyse,21:22,:,:),1),2)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\TEI_click_pow_Layer_23')
    close(fH)
    
    % Layer 4
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(data(cfg.exps2analyse,23:24,:,:),1),2)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\TEI_click_pow_Layer_24')
    close(fH)
    
    % Layer 5
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(data(cfg.exps2analyse,25:28,:,:),1),2)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\TEI_click_pow_Layer_25')
    close(fH)
    
    % Layer 6
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(data(cfg.exps2analyse,29:32,:,:),1),2)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\TEI_click_pow_Layer_26')
    close(fH)
    
    % Grand mean
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(data(cfg.exps2analyse,17:32,:,:),1),2)))
    shading interp
    caxis([-5 5])
    colormap(jet)
    print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\TEI_click_grandMean')
    close(fH)
    
    % ECoG TFR
    for ch = 33:48
        fH = figure();
        pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(data(cfg.exps2analyse,ch,:,:),1)))
        shading interp
        caxis([-5 5])
        colormap(jet)
        print(fH,'-dpng','-r1800',['C:\mbv\temp\tfr\plots\TEI_click_pow_ECoG_' num2str(ch-32)])
        close(fH)
    end
    
elseif strcmp(options.datatype,'plf')
    % Plotting phase-locking factor:
    
    % ---------
    % FIGURE 3A - Phase-locing factor - Acoustic
    % ---------
    
    % Shank 1
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(data(cfg.exps2analyse,1:16,:,:),1),2)))
    shading interp
    caxis([0 1])
    colormap(jet)
    print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\TEI_click_plf_shank1')
    close(fH)
    
    % Shank 2
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(data(cfg.exps2analyse,17:32,:,:),1),2)))
    shading interp
    caxis([0 1])
    colormap(jet)
    print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\TEI_click_plf_shank2')
    close(fH)
    
    % ECoG
    fH = figure();
    pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(data(cfg.exps2analyse,33:48,:,:),1),2)))
    shading interp
    caxis([0 1])
    colormap(jet)
    print(fH,'-dpng','-r1800','C:\mbv\temp\tfr\plots\TEI_click_plf_ECoG')
    close(fH)
end

end