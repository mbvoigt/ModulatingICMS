%% Calculate LFP - AEES = Electric only condition
%
% Subfunction for ModulatingICMS_main
%
% Calculates local field potential (LFP) responses of electrically
% stimulated trials
%
function [LFP_AEES] = Calculate_LFP_AEES(cfg)

% Initialize variable
LFP_AEES = zeros(length(cfg.experiments),length(cfg.fi),cfg.noChannels,cfg.sweepBefore+cfg.sweepAfter+1,'single');

for eID = 1:length(cfg.experiments)
    for f = 1:length(cfg.fi)
        % Loading data
        position = 1;% Added as placeholder to be able to implement different 
        % recording positions for different experiments
        pth = ['E:\ICMS-GP\Data\' cfg.experiments{eID} '\AlphaOmega\MATFiles\']; % Raw data save path
        fprintf('Single-Trial | AEES | %s | %s\n',cfg.experiments{eID},cfg.fi{f})
        expString = sprintf('Single-Trial | AEES | %s | %s | ',cfg.experiments{eID},cfg.fi{f});
        mbv_log(expString,'file',cfg.logfile);
        mbv_log([expString 'Path: ' pth],'file',cfg.logfile);
        
        % Determine the available data files to load
        if strcmp(cfg.experiments{eID},'GP14F16CL')
            folderContents=dir([pth 'ctx_pos' num2str(position) '_layering_100' cfg.fi{f} '_0001.mat']);
        else
            folderContents=dir([pth 'ctx_pos' num2str(position) '_layering_100' cfg.fi{f} '_32_0001.mat']);
        end
        mbv_log([expString 'Number of files found: ' num2str(length(folderContents))],'file',cfg.logfile);
        
        % Initialize variable
        data = zeros(cfg.noChannels,32,cfg.sweepBefore+cfg.sweepAfter+1);
        data(data==0) = NaN;
        
        % Load raw data
        load([pth folderContents(1).name],'-regexp','^CStim|^CSPK|^CTTL');
        mbv_log([expString 'File loaded: ' [pth folderContents(1).name]],'file',cfg.logfile);
        
        % Calculate trigger timestamps
        tmstmps1 = (CStimMarker_001(1,:)/(CStimMarker_001_KHz*1000))-CSPK_001_TimeBegin; %#ok Variable from file
        tm = round(tmstmps1*(CSPK_001_KHz*1000));
        mbv_log([expString 'Number of sweeps found: ' num2str(length(tm))],'file',cfg.logfile);
        
        % Cut out epoochs (=sweeps) of raw data
        for ch = 1:cfg.noChannels
            v = sprintf('CSPK_0%02d',ch);
            for sw = 1:length(tm)
                try
                    tmp = eval([v '(tm((1*' num2str(length(tm)) ')-' num2str(length(tm)) '+sw)-cfg.sweepBefore:tm((1*' num2str(length(tm)) ')-' num2str(length(tm)) '+sw)+cfg.sweepAfter)']);
                    data(ch,sw,:) = double(tmp);
                    clear tmp
                catch
                    fprintf(['Channel ' sprintf('%02d',ch) ' | ' sprintf('%02d',sw) ': Couldn''t evaluate sweep!\n'])
                    mbv_log([expString 'Channel ' num2str(ch) ' | ' num2str(sw) ': Couldn''t evaluate sweep!'],'file',cfg.logfile,'type',2);
                end
            end
        end
        
        clear('-regexp','CSPK','CStim','CTTL')
        
        % Remove first and last sweep to get rid of current source
        % artefacts (ATTN: Number of sweeps is hardcoded!)
        data = data(:,2:31,:);
        
        % Order channels according to NNX channel map
        data = data(cfg.channelidx,:,:);
        
        % Remove stimulation channel
        data(f,:,:)=NaN;
        
        % Blank stimulus artefact (ATTN: range is hardcoded here!)
        data = blankstimulus(data,22000,500,503);
        mbv_log(sprintf('%s Stimulus artefact removed by linear interpolation between 500 and 503 ms!',expString),'file',cfg.logfile);
        
        % Filter settings for local field potentials
        [z,p,k] = butter(2,150*2/22000,'low');
        [b,a] = zp2tf(z,p,k);
        
        % Initialize temporary data variable
        LFP_temp = zeros(30,cfg.sweepBefore+cfg.sweepAfter+1);
        
        for ch = 1:cfg.noChannels
            for sw = 1:30
                LFP_temp(sw,:) = single(filtfilt(b,a,squeeze(double(data(ch,sw,:))))); % 0-phase filtering
            end
            LFP_AEES(eID,f,ch,:) = squeeze(nanmean(LFP_temp));
            clear LFP_temp
        end
        
    end % repeat for every stimulated electrode
    
end % repeat for every experiment