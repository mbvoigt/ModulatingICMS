%% Calculate LFP - AEES - Pre-stimulus phase influence
%
% Subfunction for ModulatingICMS_main
%
% Calculates pre-stimulus phase influence of electrically only
% stimulated trials
%
% Plots figures 7C, 7E, 7G
%
function Calculate_LFP_AEES_PhaseInfluence(cfg)

% Initialize variables
phase = zeros(length(cfg.experiments),cfg.noChannels,30,13200,'single');
amp = zeros(length(cfg.experiments),cfg.noChannels,30);
zscore = NaN(length(cfg.experiments),cfg.noChannels,30,13200,'single');
tpoints=100:25:500; % pre-stimulus times: 500ms - tpoints 
binrange = -pi():2*pi()/8:pi(); % 8 bins, -pi to pi
MI = NaN(length(cfg.fi),length(tpoints),cfg.noChannels);

% Colormap for plots
cm = flipud([20/255 11/255 0/255;...
    100/255 55/255 0/255;...
    160/255 88/255 0/255;...
    255/255 140/255 0/255;...
    255/255 160/255 45/255;...
    255/255 205/255 145/255]);

for f = 1:length(cfg.fi)
    for eID = 1:length(cfg.experiments)
        % Load data
        position = 1;% Added as placeholder to be able to implement different 
        % recording positions for different experiments
        pth = ['E:\ICMS-GP\Data\' cfg.experiments{eID} '\AlphaOmega\MATFiles\']; % raw data save path
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
        
        % Load raw data
        load([pth folderContents(1).name],'-regexp','^CStim|^CSPK|^CTTL');
        mbv_log([expString 'File loaded: ' [pth folderContents(1).name]],'file',cfg.logfile);
        
        % Calculate trigger timestamps
        tmstmps1 = (CStimMarker_001(1,:)/(CStimMarker_001_KHz*1000))-CSPK_001_TimeBegin; %#ok Variable from file
        tm = round(tmstmps1*(CSPK_001_KHz*1000));
        mbv_log([expString 'Number of sweeps found: ' num2str(length(tm))],'file',cfg.logfile);
        
        % Cut out epoochs (=sweeps) of raw data
        data = zeros(cfg.noChannels,32,cfg.sweepBefore+cfg.sweepAfter+1); data(data==0) =NaN;
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
        
        % Get pre-stimulus data and remove first and last sweep to remove
        % edge artefacts
        data = data(:,2:31,1:cfg.sweepBefore+2200);
        
        % Reorder according to Neuronexus channel map
        data = data(cfg.channelidx,:,:);
        
        % Remove stimulation channel
        data(f,:,:)=NaN;
        
        % Blank stimulus artefact (ATTN: range is hardcoded here!)
        data = blankstimulus(data,22000,500,503);
        mbv_log(sprintf('%s Stimulus artefact removed by linear interpolation between 500 and 503 ms!',expString),'file',cfg.logfile);
        
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
                clear phasetemp
                % Z-scoring
                zscore(eID,ch,sw,:) = (LFP-(mean(LFP(1:cfg.sweepBefore))))./(std(LFP(1:cfg.sweepBefore)));
                % Determine response amplitude for each trial
                amp(eID,ch,sw) = max(abs(zscore(eID,ch,sw,cfg.sweepBefore:cfg.sweepBefore+50*22)));
            end
        end
        
        clear LFP data
        
    end % repeat for all experiments
    
    clear eID a b z p k v tm tmstmps1 sw pth position ch folderContents
    
    % Modulation index (MI) calculation
    for tidx = 1:length(tpoints)
        t = tpoints(tidx)*22;
        for ch = 1:48
            if ch ~= f
                tmpph = squeeze(phase(:,ch,:,t));
                tmpamp = squeeze(amp(:,ch,:));
                [~,~,bin] = histcounts(tmpph(:),binrange);
                tmpamp(bin==0)=[];
                bin(bin==0) = [];
                averaged = accumarray(bin, tmpamp(:),[],@mean);
                if ~isempty(averaged)
                    MI(f,tidx,ch) = range(averaged);
                end
            end % if channel is not stimulated
        end % repeat for every channel
    end % repeat for every pre-stimulus time
    
    
    % Plot single electrode examples at delta_phi = -300 ms and -5 ms:
    % Example electrode chosen for Figure 7C: electrode 7, electrode 8
    % stimulated
    for ch = 1:16
        if ch ~= f
            % Example phase influence - -5 ms
            t = 495*22;
            fH = figure('Name',['Example phase influence - -5ms - Ch. ' num2str(ch)]);
            tmpph = squeeze(phase(:,ch,:,t));
            tmpamp = squeeze(amp(:,ch,:));
            scatter(tmpph(:),tmpamp(:),'ok')
            hold on
            plot([-pi() pi()],[0 0],'--k')
            set(gca,'XTick',-pi():2*pi()/8:pi())
            [~,~,bin] = histcounts(tmpph(:),binrange);
            tmpamp(bin==0)=[];
            bin(bin==0) = [];
            averaged = accumarray(bin, tmpamp(:),[],@mean);
            averagedstd = accumarray(bin, tmpamp(:),[],@std);
            errorbar(0.5*(binrange(1:end-1)+binrange(2:end)),averaged,averagedstd,'k','LineWidth',2)
            print(fH,'-dpdf','-r1200',['C:\mbv\temp\tfr\plots\AEES_Phase\LFP_AEES_phaseScatter_-5ms_' num2str(f) '_' num2str(ch)])
            print(fH,'-dpng','-r100',['C:\mbv\temp\tfr\plots\AEES_Phase\LFP_AEES_phaseScatter_-5ms_' num2str(f) '_' num2str(ch)])
            close(fH)
            
            %
            % Example phase influence - -300 ms
            t = 200*22;
            fH = figure('Name',['Example phase influence - -300ms - Ch. ' num2str(ch)]);
            tmpph = squeeze(phase(:,ch,:,t));
            tmpamp = squeeze(amp(:,ch,:));
            scatter(tmpph(:),tmpamp(:),'ok')
            hold on
            plot([-pi() pi()],[0 0],'--k')
            set(gca,'XTick',-pi():2*pi()/8:pi())
            [~,~,bin] = histcounts(tmpph(:),binrange);
            tmpamp(bin==0)=[];
            bin(bin==0) = [];
            averaged = accumarray(bin, tmpamp(:),[],@mean);
            averagedstd = accumarray(bin, tmpamp(:),[],@std);
            errorbar(0.5*(binrange(1:end-1)+binrange(2:end)),averaged,averagedstd,'k','LineWidth',2)
            print(fH,'-dpdf','-r1200',['C:\mbv\temp\tfr\plots\AEES_Phase\LFP_AEES_phaseScatter_-300ms_' num2str(f) '_' num2str(ch)])
            print(fH,'-dpng','-r100',['C:\mbv\temp\tfr\plots\AEES_Phase\LFP_AEES_phaseScatter_-300ms_' num2str(f) '_' num2str(ch)])
            close(fH)
            
        end % plot only for electrodes which weren't stimulated
    end % repeat for every electrode
    
end % repeat for every stimulated electrode

% ---------
% FIGURE 7E - Modulation Index against delta_phi
% ---------
% Shank 1
fH = figure();
for f = 1:6
    errorbar(nanmean(nanmean(MI(cfg.layerassign{f},:,1:16)),3),nanstd(nanmean(MI(cfg.layerassign{f},:,1:16)),[],3)/sqrt(16),'Color',cm(f,:));
    hold on
end
set(gca,'YLim',[0 6])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\LFP_AEES_ModulationIndex_Time_Shank1')
close(fH)
% ---------

% ---------
% FIGURE 7G - Modulation Index against cortical layer
% ---------
% delta_phi = -300 ms
fH = figure();
for f = 1:6
    plot([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6],squeeze(nanmean(MI(cfg.layerassign{f},5,1:16))),'o','Color',cm(f,:)) % Shank 1
    hold on
    plot(accumarray([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6]', squeeze(nanmean(MI(cfg.layerassign{f},5,1:16))),[],@nanmean),'Color',cm(f,:));
end
set(gca,'XLim',[0 7],'YLim',[0 12])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\LFP_AEES_ModulationIndex_Layers_Shank1_-300ms')
close(fH)

% delta_phi = 0 ms
fH = figure();
for f = 1:6
    plot([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6],squeeze(nanmean(MI(cfg.layerassign{f},17,1:16))),'o','Color',cm(f,:)) % Shank 1
    hold on
    plot(accumarray([1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6]', squeeze(nanmean(MI(cfg.layerassign{f},17,1:16))),[],@nanmean),'Color',cm(f,:));
end
set(gca,'XLim',[0 7],'YLim',[0 12])
print(fH,'-dpdf','-r1200','C:\mbv\temp\tfr\plots\LFP_AEES_ModulationIndex_Layers_Shank1_0ms')
close(fH)
% ---------
