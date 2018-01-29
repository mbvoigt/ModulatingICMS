%% Calculate TFR - Tone-Electric Interference (TEI)
%
% Subfunction for ModulatingICMS_main
%
% Calculates time-frequency responses for auditory _and_ electrically (A&E)
% stimulated trials
%
function Calculate_TFR_TEI(cfg)

for eID = 1:length(cfg.experiments)
    
    % Loading data
    position = 1;
    % Added as placeholder to be able to implement different
    % recording positions for different experiments
    pth = ['E:\ICMS-GP\Data\' cfg.experiments{eID} '\AlphaOmega\MATFiles\']; % Raw data save path
    
    
    for d = 1 :length(cfg.delay)
        for f = 1 :length(cfg.fi)
            fprintf('TEI | %s | %s | %s\n',cfg.experiments{eID},cfg.delay{d},cfg.fi{f})
            expString = sprintf('TEI | %s | %s | %s | ',cfg.experiments{eID},cfg.delay{d},cfg.fi{f});
            mbv_log(expString,'file',cfg.logfile);
            mbv_log([expString 'Path: ' pth],'file',cfg.logfile);
            
            if ~exist(['C:\mbv\temp\tfr\TFR_' cfg.experiments{eID} '_' cfg.delay{d} '_' cfg.fi{f} '.mat'],'file') || cfg.recalculate
                
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
                data = zeros(cfg.noChannels,32,1001); 
                data(data==0) =NaN;
                
                if length(folderContents) > 1
                    
                    % Read first data file
                    load([pth folderContents(1).name],'-regexp','^CStim|^CSPK|^CTTL');
                    
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
                    load([pth folderContents(2).name],'-regexp','^CStim|^CTTL');
                   
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
                            data(ch,sw,:) = resample(double(tmp),1,22); % Resample to 1kHz!
                            clear tmp
                        end
                    end
                    
                    clear fileData
                else
                    % Load raw data
                    load([pth folderContents(1).name],'-regexp','^CStim|^CSPK|^CTTL');
                    mbv_log([expString 'File loaded: ' [pth folderContents(1).name]],'file',cfg.logfile);
                    
                    % Calculate trigger timestamps
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
                                data(ch,sw,:) = resample(double(tmp),1,22); % Resample to 1kHz!
                                clear tmp
                            catch
                                fprintf(['Channel ' sprintf('%02d',ch) ' | ' sprintf('%02d',sw) ': Couldn''t evaluate sweep!\n'])
                                mbv_log([expString 'Channel ' num2str(ch) ' | ' num2str(sw) ': Couldn''t evaluate sweep!'],'file',cfg.logfile,'type',2);
                            end
                        end
                    end
                    
                    clear('-regexp','CSPK','CStim','CTTL')
                end % handle different amount of data files
                
                % Remove first and last sweep to get rid of current source
                % artefacts (ATTN: Number of sweeps is hardcoded!)
                data = data(:,2:31,:);
                
                % Blank stimulus artefact
                data = blankstimulus(data,1000,cfg.blf(d),cfg.blt(d));
                mbv_log(sprintf('%s Stimulus artefact removed by linear interpolation between %3d and %3d ms!',expString,cfg.blf(d),cfg.blt(d)),'file',cfg.logfile);
                
                % Demean! Subtract mean of baseline from time domain data!
                for ch = 1:cfg.noChannels
                    for sw = 1:30
                        baseline = mean(squeeze(data(ch,sw,1:351)));
                        data(ch,sw,:) = data(ch,sw,:)-baseline;
                    end
                end
                mbv_log([expString 'Data demeaned by subtracting the baseline (samples 1:351) from data!'],'file',cfg.logfile);
                
                % Reformat data to fit fieldtrip data format:
                % Here data gets reformatted to a structure format,
                % at the same time channel 1:32 are now ordered according to the
                % NNX electrode order (channelidx)
                mbv_log([expString 'Reformating data to fit fieldtrip data format'],'file',cfg.logfile);
                for tr = 1:30
                    for ch=1:cfg.noChannels
                        dat.trial{tr}(ch,:) = double(data(cfg.channelidx(ch),tr,:));
                    end
                    dat.time{tr} = -0.5:1/1000:0.5;
                end
                for ch = 1:cfg.noChannels
                    dat.label{ch}=num2str(ch);
                end
                dat.fsample=1000;
                clear tr ch data
                
                % Log FT configuration
                cfgnames = fieldnames(cfg.FTcfg);
                for idx = 1:length(cfgnames)
                    if strcmp('foi',cfgnames{idx})
                        mbv_log([expString 'FieldTrip configuration: foi : 7:2:95'],'file',cfg.logfile);
                    elseif strcmp('toi',cfgnames{idx})
                        mbv_log([expString 'FieldTrip configuration: toi : -0.5:0.001:0.5'],'file',cfg.logfile);
                    else
                        mbv_log([expString 'FieldTrip configuration: ' sprintf('%s : %s',cfgnames{idx},num2str(cfg.FTcfg.(cfgnames{idx})))],'file',cfg.logfile);
                    end
                end
                clear idx cfgnames
                
                % Time-Frequency Analysis
                mbv_log([expString 'Calculating time frequency representation using fieldtrip'],'file',cfg.logfile);
                TFR = ft_freqanalysis(cfg.FTcfg,dat);
                fprintf(repmat('\b',1,826)) % Remove output of freqanalysis
                fprintf('\n')
                clear dat
                
                % Extracting the complex fourier spectrum
                fspctrm = single(TFR.fourierspctrm);
                TFR = rmfield(TFR,'fourierspctrm');
                
                % phase-locking factor
                plfspctrm = fspctrm./abs(fspctrm);
                plfspctrm = squeeze(nanmean(plfspctrm));
                plfspctrm = abs(plfspctrm);
                
                % power spectrum
                powspctrm = abs(fspctrm).^2;
                clear fspctrm
                
                % Single-Trial-Low-Gamma-Power (ST_LGP) calculation for influence
                % of pre-stimulus phase
                clear temp ST_LGP
                temp = zeros(length(cfg.FTcfg.foi),1001);
                ST_LGP = zeros(cfg.noChannels,30);
                ST_LGPE = zeros(cfg.noChannels,30);
                for sw = 1:30
                    for ch = 1:cfg.noChannels
                        % Baseline subtraction
                        for freq = 1:length(cfg.FTcfg.foi)
                            baseline = squeeze(nanmean(powspctrm(sw,ch,freq,cfg.BaselineWindow)));
                            temp(freq,:) = 10*log10(powspctrm(sw,ch,freq,:)/baseline);
                        end
                        % Extraction of low-gamma power from analysis
                        % window
                        ST_LGP(ch,sw) = max(nanmean(temp(cfg.FTcfg.foi>29&cfg.FTcfg.foi<60,cfg.InducedResponse),1),[],2);
                        ST_LGPE(ch,sw) = max(nanmean(temp(cfg.FTcfg.foi>29&cfg.FTcfg.foi<60,cfg.EvokedResponse),1),[],2);
                    end
                end
                mbv_log('Extracted single-trial low-gammma power values (ST_LGP)','file',cfg.logfile);
                %
                
                % Take the _median_ over single trials
                powspctrm = squeeze(nanmedian(powspctrm));
                
                % Save processed data to file
                TFR.dimord = TFR.dimord(5:end);
                TFR.plfspctrm = plfspctrm; clear plfspctrm
                TFR.powspctrm = powspctrm; clear powspctrm
                save(['C:\mbv\temp\tfr\TFR_' cfg.experiments{eID} '_' cfg.delay{d} '_' cfg.fi{f}],'TFR','ST_LGP','ST_LGPE')
                mbv_log([expString 'Results saved at: ' 'C:\mbv\temp\tfr\TFR_' cfg.experiments{eID} '_' cfg.delay{d} '_' cfg.fi{f} '.mat'],'file',cfg.logfile);
                
                clear TFR
                
            end % if processed data file already exists and no recalculation is whished for, skip everything!
        end % repeat for every stimulated electrode
    end % repeat for every stimulus delay
end % repeat for every experiment
