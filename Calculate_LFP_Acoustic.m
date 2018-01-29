%% Calculate LFP - Acoustic
%
% Subfunction for ModulatingICMS_main
%
% Calculates local field potential (LFP) responses of acoustically 
% stimulated trials
%
function LFP = Calculate_LFP_Acoustic(cfg)

% Initialize variable
LFP = zeros(length(cfg.experiments),cfg.noChannels,30,cfg.sweepBefore+cfg.sweepAfter+1);

for eID = 1:length(cfg.experiments)
    % Loading data
    position = 1; % Added as placeholder to be able to implement different 
    % recording positions for different experiments
    pth = ['E:\ICMS-GP\Data\' cfg.experiments{eID} '\AlphaOmega\MATFiles\']; % Raw data save path
    fprintf('LFP | %s | Acoustic | \n',cfg.experiments{eID})
    expString = sprintf('LFP | %s | Acoustic | ',cfg.experiments{eID});
    mbv_log(expString,'file',cfg.logfile);
    mbv_log([expString 'Path: ' pth],'file',cfg.logfile);
    
    % Determine the available data files to load
    folderContents=dir([pth 'ctx_pos' num2str(position) '_click00' cfg.clickfile{eID} '.mat']);
    mbv_log([expString 'Number of files found: ' num2str(length(folderContents))],'file',cfg.logfile);
    
    % Initialize data variable
    data = zeros(cfg.noChannels,30,22001);
    data(data==0) = NaN;
    
    % Load raw data
    load([pth folderContents(1).name],'-regexp','^CStim|^CSPK_|^CTTL');
    mbv_log([expString 'File loaded: ' [pth folderContents(1).name]],'file',cfg.logfile);
    
    % Calculate trigger timestamps
    timeoffset = (CTTL_001_TimeBegin)-(CSPK_001_TimeBegin);
    tmstmps1 = (CTTL_001_Down(1,:)/(CTTL_001_KHz*1000))+timeoffset; %#ok Variable is loaded from file
    tm = round(tmstmps1*(CSPK_001_KHz*1000));
    clear tmstmps1 timeoffset
    mbv_log([expString 'Number of sweeps found in file: ' num2str(length(tm))],'file',cfg.logfile);
    
    % Cut out epoochs (=sweeps) of raw data
    for ch = 1:cfg.noChannels
        v =  sprintf('CSPK_0%02d',ch);
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
    
    % Filter settings for local field potentials
    [z,p,k] = butter(2,150*2/22000,'low');
    [b,a] = zp2tf(z,p,k);
    
    % Initialize filtered data variable
    dataFiltered = zeros(cfg.noChannels,30,22001);
    for ch = 1:cfg.noChannels
        for sw = 1:30
            dataFiltered(ch,sw,:) = filtfilt(b,a,data(ch,sw,:)); % 0-phase filtering
        end
    end
    clear data
    

    % Order channels according to NNX channel map
    LFP(eID,:,:,:) = dataFiltered(cfg.channelidx,:,:);
    clear dataFiltered
    
end % repeat for every experiment

% Clear current_monitor channel no 33 from experiments 1 & 2
LFP([1,2],33,:,:) = NaN;
