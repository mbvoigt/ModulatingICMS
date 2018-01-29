%% Calculate LFP - TEI - Pre-stimulus phase influence
%
% Subfunction for ModulatingICMS_main
%
% Calculates pre-stimulus phase influence of acoustically _and_ electrically
% stimulated trials
%
% Plots figures 7H, 7I, 7J
%
function Calculate_LFP_TEI_PhaseInfluence(cfg)

% Initialize variables
phase = zeros(length(cfg.experiments),cfg.noChannels,30,13200,'single');
tpoints=100:25:500; % pre-stimulus times: 500ms - tpoints
binrange = -pi():2*pi()/8:pi(); % 8 bins, -pi to pi
amp = zeros(length(cfg.experiments),cfg.noChannels,30);
ampE = zeros(length(cfg.experiments),cfg.noChannels,30);
MI = NaN(length(cfg.delay),length(cfg.fi),length(tpoints),cfg.noChannels); % Induced response
MIE = NaN(length(cfg.delay),length(cfg.fi),length(tpoints),cfg.noChannels); % Evoked resposne

for d = 1:length(cfg.delay)
    for f = 1:length(cfg.fi)
        for eID = 1:length(cfg.experiments)
            % Load data
            position = 1;% Added as placeholder to be able to implement different
            % recording positions for different experiments
            pth = ['E:\ICMS-GP\Data\' cfg.experiments{eID} '\AlphaOmega\MATFiles\']; % raw data save path
            fprintf('TEI | %s | %s | %s\n',cfg.experiments{eID},cfg.delay{d},cfg.fi{f})
            expString = sprintf('TEI | %s | %s | %s | ',cfg.experiments{eID},cfg.delay{d},cfg.fi{f});
            fprintf([expString '\n'])
            mbv_log(expString,'file',cfg.logfile);
            mbv_log([expString 'Path: ' pth],'file',cfg.logfile);
            
            % Determine the available data files to load
            if strcmp(cfg.experiments{eID},'GP17K16CL') && strcmp(cfg.delay(d),'5')
                folderContents=dir([pth 'ctx_pos' num2str(position) '_tei_0' cfg.delay{d} 'ms_100' cfg.fi{f} '0001.mat']);
            elseif strcmp(cfg.experiments{eID},'GP22K16CL') && strcmp(cfg.delay(d),'5')
                folderContents=dir([pth 'ctx_pos' num2str(position) '_tei_0' cfg.delay{d} 'ms_100' cfg.fi{f} '0001.mat']);
            elseif strcmp(cfg.experiments{eID},'GP24K16CL') && strcmp(cfg.delay(d),'5')
                folderContents=dir([pth 'ctx_pos' num2str(position) '_tei_0' cfg.delay{d} 'ms_100' cfg.fi{f} '0001.mat']);
            elseif strcmp(cfg.experiments{eID},'GP17K16CL') && strcmp(cfg.delay(d),'15')
                folderContents=dir([pth 'ctx_pos' num2str(position) '_tei_' cfg.delay{d} '_ms_100' cfg.fi{f} '0001.mat']);
            else
                folderContents=dir([pth 'ctx_pos' num2str(position) '_tei_' cfg.delay{d} 'ms_100' cfg.fi{f} '0001.mat']);
            end
            mbv_log([expString 'Number of files found: ' num2str(length(folderContents))],'file',cfg.logfile);
            
            % Initialize data variable
            data = zeros(cfg.noChannels,32,cfg.sweepBefore+cfg.sweepAfter+1);
            data(data==0) =NaN;
            
            if length(folderContents) > 1
                % Read first data file
                load([pth folderContents(1).name],'-regexp','^CStim|^CSPK_00|^CSPK_01|^CTTL');
                % Process timestamps of first file
                timeoffset = (CTTL_001_TimeBegin)-(CSPK_001_TimeBegin);
                tmstmps1 = (CTTL_001_Down(1,:)/(CTTL_001_KHz*1000))+timeoffset; %#ok variable is loaded from file
                tm1 = round(tmstmps1*(CSPK_001_KHz*1000));
                offset = length(CSPK_001);
                % Prepare raw data variable for first file
                for ch = 1:cfg.noChannels
                    v =  sprintf('CSPK_0%02d',ch);
                    fileData.(v) = int16(eval(v));
                end
                clear('-regexp','^CSPK|^CStim|^CTTL');
                % Read second data file
                load([pth folderContents(2).name],'-regexp','^CStim|^CSPK_00|^CSPK_01|^CTTL');
                % Process timestamps of second file
                timeoffset = (CTTL_001_TimeBegin)-(CSPK_001_TimeBegin);
                tmstmps2 = (CTTL_001_Down(1,:)/(CTTL_001_KHz*1000))+timeoffset;
                tm2 = round(tmstmps2*(CSPK_001_KHz*1000));
                % Add raw data from second file to prepared raw data variable
                for ch = 1:cfg.noChannels
                    v = sprintf('CSPK_0%02d',ch);
                    fileData.(v) = [fileData.(v) int16(eval(v))];
                end
                clear('-regexp','^CSPK|^CStim|^CTTL');
                % Combine timestamps of both files
                if exist('tm2','var') == 1
                    tm = [tm1 tm2+offset];
                else
                    tm = tm1;
                end
                clear tm1 tm2 tmstmps1 tmstmps2 timeoffset v ch
                
                % Cut out epoochs (=sweeps) of raw data
                for ch = 1:cfg.noChannels
                    v = sprintf('CSPK_0%02d',ch); %#ok Value used by eval
                    for sw = 1:length(tm)
                        tmp = eval(['fileData.(v)(tm((1*' num2str(length(tm)) ')-' num2str(length(tm)) '+sw)-cfg.sweepBefore:tm((1*' num2str(length(tm)) ')-' num2str(length(tm)) '+sw)+cfg.sweepAfter)']);
                        data(ch,sw,:) = double(tmp);
                        clear tmp
                    end
                end
                
                clear fileData
            else % if only one file to load
                % Read data file
                load([pth folderContents(1).name],'-regexp','^CStim|^CSPK_00|^CSPK_01|^CSPK_02|^CSPK_03|^CTTL');
                mbv_log([expString 'File loaded: ' [pth folderContents(1).name]],'file',cfg.logfile);
                
                % Process timestamps
                timeoffset = (CTTL_001_TimeBegin)-(CSPK_001_TimeBegin);
                tmstmps1 = (CTTL_001_Down(1,:)/(CTTL_001_KHz*1000))+timeoffset; %#ok variable is loaded from file
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
            end
            
            % Remove first and last sweep to remove edge artefacts
            data = data(:,2:31,:);
            
            % Order channels according to channel map
            data = data(cfg.channelidx,:,:);
            
            % Remove stimulation channel
            data(f,:,:) = NaN;
            
            % Blank stimulus artefact
            data = blankstimulus(data,22000,cfg.blf(d),cfg.blt(d));
            mbv_log(sprintf('%s Stimulus artefact removed by linear interpolation between %3d and %3d ms!',expString,cfg.blf(d),cfg.blt(d)),'file',cfg.logfile);
            
            % Filter settings for local field potentials
            [z,p,k] = butter(2,150*2/22000,'low');
            [b,a] = zp2tf(z,p,k);
            
            % Filtering and phase calculation
            for ch = 1:cfg.noChannels
                for sw = 1:30
                    % LFP filtering
                    LFP = filtfilt(b,a,squeeze(double(data(ch,sw,:))));
                    % Phase extraction
                    phasetemp = hilbert(LFP);
                    phase(eID,ch,sw,:) = angle(phasetemp);
                end
            end
            clear data phasetemp LFP
            
            % Load corresponding single-trial TFR amplitude data
            load(['C:\mbv\temp\tfr\TFR_' cfg.experiments{eID} '_' cfg.delay{d} '_' cfg.fi{f} '.mat'],'ST_LGP')
            % Induced response
            amp(eID,:,:) = ST_LGP;
            % Evoked response
            load(['C:\mbv\temp\tfr\TFR_' cfg.experiments{eID} '_' cfg.delay{d} '_' cfg.fi{f} '.mat'],'ST_LGPE')
            ampE(eID,:,:) = ST_LGPE;
            clear ST_LGP ST_LGPE
            
        end % repeat for all experiments
        
        fprintf([expString ' Calculating phase influence\n'])
        
        % Modulation index calculation
        for tidx = 1:length(tpoints)
            t = tpoints(tidx)*22;
            for ch = 1:48
                if ch ~= f
                    tmpph = squeeze(phase(:,ch,:,t));
                    tmpamp = squeeze(amp(:,ch,:));
                    tmpampE = squeeze(ampE(:,ch,:));
                    [~,~,bin] = histcounts(tmpph(:),binrange);
                    tmpamp(bin==0)=[];
                    tmpampE(bin==0)=[];
                    bin(bin==0) = [];
                    % Induced response
                    averaged = accumarray(bin, tmpamp(:),[],@mean);
                    if ~isempty(averaged)
                        MI(d,f,tidx,ch) = range(averaged);
                    end
                    % Evoked response
                    averagedE = accumarray(bin, tmpampE(:),[],@mean);
                    if ~isempty(averagedE)
                        MIE(d,f,tidx,ch) = range(averagedE);
                    end
                    
                end % if channel is not stimulated
            end % repeat for every channel
        end % repeat for every pre-stimulus time
        
    end % repeat for every electrode stimulated
end % repeat for every stimulus delay

% Save processing to file
save('C:\mbv\temp\tfr\STP_TEI','MI','Explained','tpoints','MIE','ExplainedE')

%%
% Load processing from file
% load('C:\mbv\temp\tfr\STP_TEI','MI','Explained','tpoints','MIE','ExplainedE')

% Colormap for plots
cm = flipud([20/255 11/255 0/255;...
    100/255 55/255 0/255;...
    160/255 88/255 0/255;...
    255/255 140/255 0/255;...
    255/255 160/255 45/255;...
    255/255 205/255 145/255]);

for d=1:3 % for each stimulus delay
    % Stimulus delay used for Figure 7: d = 1 -> 5 ms
    %
    % ---------
    % FIGURE 7I - Modulation Index against delta_phi, induced response
    % ---------
    % Shank 1
    fH = figure();
    for f = 1:6
        errorbar(squeeze(nanmean(nanmean(MI(d,cfg.layerassign{f},9:17,1:16),2),4)),squeeze(nanstd(nanmean(MI(d,cfg.layerassign{f},9:17,1:16),2),[],4)/sqrt(16)),'Color',cm(f,:));
        hold on
    end
    set(gca,'YLim',[0 7])
    print(fH,'-dpdf','-r1200',['C:\mbv\temp\tfr\plots\TEI_STP_ModulationIndex_Time_Shank1_' num2str(d)])
    close(fH)
    
    % ---------
    % FIGURE 7H - Modulation Index against cortical layer, evoked response
    % ---------
    % delta_phi = 0 ms
    fH = figure();
    for f = 1:6
        plot([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6],squeeze(nanmean(MIE(d,cfg.layerassign{f},17,1:16),2)),'o','Color',cm(f,:)) % Shank 1
        hold on
        plot(accumarray([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6]', squeeze(nanmean(MIE(d,cfg.layerassign{f},17,1:16),2)),[],@nanmean),'Color',cm(f,:));
    end
    set(gca,'XLim',[0 7],'YLim',[0 18])
    print(fH,'-dpdf','-r1200',['C:\mbv\temp\tfr\plots\TEI_STP_Evoked_ModulationIndex_Layers_Shank1_0ms_' num2str(d)])
    close(fH)
    
    % delta_phi = -300 ms
    fH = figure();
    for f = 1:6
        plot([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6],squeeze(nanmean(MIE(d,cfg.layerassign{f},5,1:16),2)),'o','Color',cm(f,:)) % Shank 1
        hold on
        plot(accumarray([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6]', squeeze(nanmean(MIE(d,cfg.layerassign{f},5,1:16),2)),[],@nanmean),'Color',cm(f,:));
    end
    set(gca,'XLim',[0 7],'YLim',[0 18])
    print(fH,'-dpdf','-r1200',['C:\mbv\temp\tfr\plots\TEI_STP_Evoked_ModulationIndex_Layers_Shank1_-300ms_' num2str(d)])
    close(fH)
    
    % ---------
    % FIGURE 7J - Modulation Index against cortical layer, induced response
    % ---------
    % delta_phi = 0 ms
    fH = figure();
    for f = 1:6
        plot([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6],squeeze(nanmean(MI(d,cfg.layerassign{f},17,1:16),2)),'o','Color',cm(f,:)) % Shank 1
        hold on
        plot(accumarray([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6]', squeeze(nanmean(MI(d,cfg.layerassign{f},17,1:16),2)),[],@nanmean),'Color',cm(f,:));
    end
    set(gca,'XLim',[0 8],'YLim',[0 15])
    print(fH,'-dpdf','-r1200',['C:\mbv\temp\tfr\plots\TEI_STP_Induced_ModulationIndex_Layers_Shank1_0ms_' num2str(d)])
    close(fH)
    
    % delta_phi = -300 ms
    fH = figure();
    for f = 1:6
        plot([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6],squeeze(nanmean(MI(d,cfg.layerassign{f},5,1:16),2)),'o','Color',cm(f,:)) % Shank 1
        hold on
        plot(accumarray([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6]', squeeze(nanmean(MI(d,cfg.layerassign{f},5,1:16),2)),[],@nanmean),'Color',cm(f,:));
    end
    set(gca,'XLim',[0 8],'YLim',[0 15])
    print(fH,'-dpdf','-r1200',['C:\mbv\temp\tfr\plots\TEI_STP_Induced_ModulationIndex_Layers_Shank1_-300ms_' num2str(d)])
    close(fH)
    
end % repeat for each stimulus delay

