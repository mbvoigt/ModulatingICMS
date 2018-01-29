%% Load TFR - TEI (= Tone-Electric Interference, combined stimulation)
%
% Subfunction for ModulatingICMS_main
%
% Loads TFRs for combined acosutically _and_ electrically stimulated trials
%
function varargout = Load_TFR_TEI(cfg,varargin)
% call with parameter 'save2file' set to 1 to save processed TFRs to disk
% called without parameter 'save2file' (or set to 0, default) loads
% processed TFRs to output variables

% Check function input
if length(varargin) == 2
    pair = reshape(varargin,2,[]);
    if strcmp(pair{1},'save2file')
        save2file = pair{2};
    else
        disp('Only save2file and a logical value is supported as name-value pair!')
        return
    end
elseif length(varargin) ~= 2 && ~isempty(varargin)
    disp('Only one name-value pair consisting of save2file and a logical value is supported!')
    return
else
    save2file = 0;
end

% -----------------------------
% Load & process delta_t = 5 ms
% -----------------------------

% Initialize variables
pow_5ms(length(cfg.experiments),length(cfg.fi),length(cfg.ch2analyse),length(cfg.FTcfg.foi),1001) = 0;
nRemoved5 = 0; % count noisy channels

for eID = 1:length(cfg.experiments)
    for f = 1:length(cfg.fi)
        
        fprintf('TEI-Analysis| %s | TEI | 5ms| Stim. electrode %s | Loading TFR\n',cfg.experiments{eID},cfg.fi{f})
        mbv_log(sprintf('TEI-Analysis| %s | TEI | 5ms | Stim. electrode %s | Loading TFR',cfg.experiments{eID},cfg.fi{f}),'file',cfg.logfile);
        % Load TFR data
        load(['C:\mbv\temp\tfr\TFR_' cfg.experiments{eID} '_5_' cfg.fi{f}]);
        
        % Get power spectrum
        pow_5ms(eID,f,:,:,:) = TFR.powspctrm(:,:,:);
        clear TFR
        
        % Remove stimulation channel
        pow_5ms(eID,f,f,:,:) = NaN;
        
        % Baseline subtraction
        for ch = 1:length(cfg.ch2analyse)
            for freq = 1:length(cfg.FTcfg.foi)
                baseline = squeeze(nanmean(pow_5ms(eID,f,ch,freq,cfg.BaselineWindow)));
                pow_5ms(eID,f,ch,freq,:) = 10*log10(pow_5ms(eID,f,ch,freq,:)/baseline);
            end
        end
        
        % Remove channels with too much line noise = 50 Hz
        for ch = 1:length(cfg.ch2analyse)
            if nanmean(nanmean(pow_5ms(eID,f,ch,cfg.FTcfg.foi>45&cfg.FTcfg.foi<55,500:end)))<-6
                pow_5ms(eID,f,ch,:,:) = NaN;
                nRemoved5 = nRemoved5+1;
            end
        end
        
    end % repeat for every stimulated electrode
end % repeat for every experiment

% Remove channel 33 from experiments 1&2
% (Channel 33 was the current monitor resistor in experiments 1&2)
try
    pow_5ms([1,2],:,33,:,:) = NaN;
catch
    fprintf('Tried to remove channel 33 from experiments 1&2: Unsuccessful!\n')
end

% Either save to file or pass on as output
if save2file
    save('C:\mbv\temp\TEI_analyses_TEI_5ms_pow','pow_5ms','nRemoved5')
else
    varargout{1} = pow_5ms;
end
clear pow_5ms nRemoved5

% ------------------------------
% Load & process delta_t = 15 ms
% ------------------------------

% Initalize variables
pow_15ms(length(cfg.experiments),length(cfg.fi),length(cfg.ch2analyse),length(cfg.FTcfg.foi),1001) = 0;
nRemoved15 = 0; % count noisy channels

for eID = 1:length(cfg.experiments)
    for f = 1:length(cfg.fi)
        
        fprintf('TEI-Analysis| %s | TEI | 15ms| Stim. electrode %s | Loading TFR\n',cfg.experiments{eID},cfg.fi{f})
        mbv_log(sprintf('TEI-Analysis| %s | TEI | 15ms | Stim. electrode %s | Loading TFR',cfg.experiments{eID},cfg.fi{f}),'file',cfg.logfile);
        % Load TFR data
        load(['C:\mbv\temp\tfr\TFR_' cfg.experiments{eID} '_15_' cfg.fi{f}]);
        
        % Get power spectrum
        pow_15ms(eID,f,:,:,:) = TFR.powspctrm(:,:,:);
        clear TFR
        
        % Remove stimulation channel
        pow_15ms(eID,f,f,:,:) = NaN;
        
        
        % Baseline subtraction
        for ch = 1:length(cfg.ch2analyse)
            for freq = 1:length(cfg.FTcfg.foi)
                baseline = squeeze(nanmean(pow_15ms(eID,f,ch,freq,cfg.BaselineWindow)));
                pow_15ms(eID,f,ch,freq,:) = 10*log10(pow_15ms(eID,f,ch,freq,:)/baseline);
            end
        end
        
        % Remove channels with too much line noise = 50 Hz
        for ch = 1:length(cfg.ch2analyse)
            if nanmean(nanmean(pow_15ms(eID,f,ch,cfg.FTcfg.foi>45&cfg.FTcfg.foi<55,500:end)))<-6
                pow_15ms(eID,f,ch,:,:) = NaN;
                nRemoved15 = nRemoved15+1;
            end
        end
        
    end % repeat for every stimulated electrode
end % repeat for every experiment

% Remove channel 33 from experiments 1&2
% (Channel 33 was the current monitor resistor in experiments 1&2)
try
    pow_15ms([1,2],:,33,:,:) = NaN;
catch
    fprintf('Tried to remove channel 33 from experiments 1&2: Unsuccessful!\n')
end

% Either save to file or pass on as output
if save2file
    save('C:\mbv\temp\TEI_analyses_TEI_15ms_pow','pow_15ms','nRemoved15')
else
    varargout{2} = pow_15ms;
end
clear pow_15ms nRemoved15

% ------------------------------
% Load & process delta_t = 25 ms
% ------------------------------

% Initialize variables
pow_25ms(length(cfg.experiments),length(cfg.fi),length(cfg.ch2analyse),length(cfg.FTcfg.foi),1001) = 0;
nRemoved25 = 0;% count noisy channels

for eID = 1:length(cfg.experiments)
    for f = 1:length(cfg.fi)
        
        fprintf('TEI-Analysis| %s | TEI | 25ms| Stim. electrode %s | Loading TFR\n',cfg.experiments{eID},cfg.fi{f})
        mbv_log(sprintf('TEI-Analysis| %s | TEI | 25ms | Stim. electrode %s | Loading TFR',cfg.experiments{eID},cfg.fi{f}),'file',cfg.logfile);
        % Load TFR data
        load(['C:\mbv\temp\tfr\TFR_' cfg.experiments{eID} '_25_' cfg.fi{f}]);
        
        % Get power spectrum
        pow_25ms(eID,f,cfg.ch2analyse,:,:) = TFR.powspctrm(cfg.ch2analyse,:,:);
        clear TFR
        
        % Remove stimulation channel
        pow_25ms(eID,f,f,:,:) = NaN;
        
        % Baseline subtraction
        for ch = 1:length(cfg.ch2analyse)
            for freq = 1:length(cfg.FTcfg.foi)
                baseline = squeeze(nanmean(pow_25ms(eID,f,ch,freq,cfg.BaselineWindow)));
                pow_25ms(eID,f,ch,freq,:) = 10*log10(pow_25ms(eID,f,ch,freq,:)/baseline);
            end
        end
        
        % Remove channels with too much line noise = 50 Hz
        for ch = 1:length(cfg.ch2analyse)
            if nanmean(nanmean(pow_25ms(eID,f,ch,cfg.FTcfg.foi>45&cfg.FTcfg.foi<55,500:end)))<-6
                pow_25ms(eID,f,ch,:,:) = NaN;
                nRemoved25 = nRemoved25+1;
            end
        end
        
    end % repeat for every stimulated electrode
end % repeat for every experiment

% Remove channel 33 from experiments 1&2
% (Channel 33 was the current monitor resistor in experiments 1&2)
try
    pow_25ms([1,2],:,33,:,:) = NaN;
catch
    fprintf('Tried to remove channel 33 from experiments 1&2: Unsuccessful!\n')
end

% Either save to file or pass on as output
if save2file
    save('C:\mbv\temp\TEI_analyses_TEI_25ms_pow','pow_25ms','nRemoved25')
else
    varargout{3} = pow_25ms;
end