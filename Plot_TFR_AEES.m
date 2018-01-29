%% Plot TFR - AEES
%
% Subfunction for ModulatingICMS_main
%
% Plots TFRs for electrically stimulated trials
%
% Plots figures 4A, 4B
%
function Plot_TFR_AEES(cfg,data,varargin)
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
    
    % Reduce 16 stimElectrodes to 6 stimulated layers
    temp(1,:,:,:) = squeeze(nanmean(nanmean(data(cfg.exps2analyse,1:2,cfg.ch2analyse,:,:),2),1));
    temp(2,:,:,:) = squeeze(nanmean(nanmean(data(cfg.exps2analyse,3:4,cfg.ch2analyse,:,:),2),1));
    temp(3,:,:,:) = squeeze(nanmean(nanmean(data(cfg.exps2analyse,5:6,cfg.ch2analyse,:,:),2),1));
    temp(4,:,:,:) = squeeze(nanmean(nanmean(data(cfg.exps2analyse,7:8,cfg.ch2analyse,:,:),2),1));
    temp(5,:,:,:) = squeeze(nanmean(nanmean(data(cfg.exps2analyse,9:12,cfg.ch2analyse,:,:),2),1));
    temp(6,:,:,:) = squeeze(nanmean(nanmean(data(cfg.exps2analyse,13:16,cfg.ch2analyse,:,:),2),1));
    
    % ---------
    % FIGURE 4A - Time-frequency response
    % ---------
    for stimL = 1:6
        % Shank 1
        fH = figure();
        pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(temp(stimL,1:16,:,:),2)))
        shading interp
        caxis([-5 5])
        colormap(jet)
        print(fH,'-dpng','-r1800',['C:\mbv\temp\tfr\plots\TEI_AEES_pow_Layer_1' num2str(stimL)])
        close(fH)
        
        % Shank 2
        fH = figure();
        pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(temp(stimL,17:32,:,:),2)))
        shading interp
        caxis([-5 5])
        colormap(jet)
        print(fH,'-dpng','-r1800',['C:\mbv\temp\tfr\plots\TEI_AEES_pow_Layer_2' num2str(stimL)])
    end
    % ---------

elseif strcmp(options.datatype,'plf')
    % Plotting phase-locking factor:
    
    % Reduce 16 stimElectrodes to 6 stimulated layers
    temp(1,:,:,:) = squeeze(nanmean(nanmean(data(cfg.exps2analyse,1:2,cfg.ch2analyse,:,:),2),1));
    temp(2,:,:,:) = squeeze(nanmean(nanmean(data(cfg.exps2analyse,3:4,cfg.ch2analyse,:,:),2),1));
    temp(3,:,:,:) = squeeze(nanmean(nanmean(data(cfg.exps2analyse,5:6,cfg.ch2analyse,:,:),2),1));
    temp(4,:,:,:) = squeeze(nanmean(nanmean(data(cfg.exps2analyse,7:8,cfg.ch2analyse,:,:),2),1));
    temp(5,:,:,:) = squeeze(nanmean(nanmean(data(cfg.exps2analyse,9:12,cfg.ch2analyse,:,:),2),1));
    temp(6,:,:,:) = squeeze(nanmean(nanmean(data(cfg.exps2analyse,13:16,cfg.ch2analyse,:,:),2),1));
    
    % ---------
    % FIGURE 4B - Phase-locking factor
    % ---------
    for stimL = 1:6
        % Shank 1
        fH = figure();
        pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(temp(stimL,1:16,:,:),1),2)))
        shading interp
        caxis([0 1])
        colormap(jet)
        print(fH,'-dpng','-r1800',['C:\mbv\temp\tfr\plots\TEI_AEES_plf_shank1_stimL_' num2str(stimL)])
        close(fH)
        
        % Shank 2
        fH = figure();
        pcolor(-0.5:0.001:0.5,7:2:95,squeeze(nanmean(nanmean(temp(stimL,17:32,:,:),1),2)))
        shading interp
        caxis([0 1])
        colormap(jet)
        print(fH,'-dpng','-r1800',['C:\mbv\temp\tfr\plots\TEI_AEES_plf_shank2_stimL_' num2str(stimL)])
        close(fH)
    end
    % ---------
    
end
