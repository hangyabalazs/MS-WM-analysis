function [spsth_all, stats_all, time, spt] = ultimate_psth_wm(cellids, event_type, event,window,varargin)
%ULTIMATE_PSTH   Peri-stimulus time histogram.
%   [PSTH SPSTH SPSTH_SE] = ULTIMATE_PSTH(CELLID,EVENT_TYPE,EVENT,WINDOW,VARARGIN)
%   calculates peri-stimulus time histogram (PSTH) for the cell passed in
%   CELLID. Smoothed PSTH (SPSTH) and SE of smoothed PSTH (SPSTH_SE) are
%   also returned.
%
%   [PSTH SPSTH SPSTH_SE TAGS] = ULTIMATE_PSTH(CELLID,EVENT_TYPE,EVENT,WINDOW,VARARGIN)
%   returns partition tags (TAGS) corrsponding to PSTHs when trials are
%   partitioned; see PARTITION_TRIALS.
%
%   [PSTH SPSTH SPSTH_SE TAGS SPT] = ULTIMATE_PSTH(CELLID,EVENT_TYPE,EVENT,WINDOW,VARARGIN)
%   returns the bin raster (SPT); see STIMES2BINRASTER.
%
%   [PSTH SPSTH SPSTH_SE SPT TAGS STATS] = ULTIMATE_PSTH(CELLID,EVENT_TYPE,EVENT,WINDOW,VARARGIN)
%   calculates and returns test results for significant firing rate changes
%   after the event (see PSTH_STATS for details).
%
%   ULTIMATE_PSTH is also capable of using two different events for the
%   periods before and after 0, usefull for statistical testing with a
%   baseline period aligned to a different event than the test period (see
%   below and PSTH_STATS).
%
%   Mandatory input arguments:
%       CELLID: defines the cell (see CellBase documentation) or session
%           (for lick-PSTH)
%       EVENT: the event to which the PSTH is aligned; if EVENT is a cell
%           array of two strings, the first event is used for the PSTH
%           and binraster before 0 and the second event is used for the
%           PSTH and binraster after 0; if EVENT is a function handle, the
%           function is called for CELLID to define the aligning event
%           (dynamic event definition)
%       EVENT_TYPE: the type of event, 'trial' 
%       WINDOW: window for calculation relative to the event in seconds
%
%   Default behavior of ULTIMATE_PSTH can be modified by using a set of
%   paramter-value pairs as optional input parameters. The following
%   parameters are implemented (with default values):
%   	'dt', 0.001 - time resolution in seconds
%       'sigma', 0.02 - smoothing kernel for the smoothed PSTH, in seconds
%       'margin',[-0.01 0.01] margins for PSTH calculation to get rid of
%           edge effect due to smoothing
%       'event_filter', 'none' - filter light-stimulation trials; see
%           FILTERTRIALS for implemented filter types
%       'filterinput',[] - some filters require additional input; see
%           FILTERTRIALS for details
%       'maxtrialno', 5000 - maximal number of trials included; if ther are
%           more valid trials, they are randomly down-sampled
%       'data_type', 'real' - ultimate_psth handles virtual/simulated data;
%           set this property to 'virtual' to analyze simulated spikes
%       'first_event', [] - event name used to exclude spikes before
%           previous event
%       'last_event', [] - event name used to exclude spikes after
%           following event
%       'spike_def' [] - if it is set to 'burst', only bursts, when set to
%           'single', only single spikes are used
%       'parts', 'all' - partitioning the set of trials; input to
%           PARTITION_TRIALS, see details therein (default, no
%           partitioning)
%       'isadaptive', 1 - 0, classic PSTH algorithm is applied; 1, adaptive
%           PSTH is calculated (see APSTH); 2, 'doubly adaptive' PSTH
%           algorithm is used (see DAPSTH)
%       'forcesmoothedstat', false - if true, smoothe PSTH is used for
%           calculating PSTH_STATS, even in the adaptive cases
%   	'baselinewin', [-0.25 0] - limits of baseline window for
%           statistical testing (see PSTH_STATS), time relative to 0 in
%           seconds
%   	'testwin', [0 0.1] - limits of test window for statistical testing
%           (see PSTH_STATS), time relative to 0 in seconds
%       'relative_threshold', 0.5 - threshold used to assess start and end
%           points of activation and inhibition intervals in PSTH_STATS; in
%           proportion of the peak-baseline difference (see PSTH_STATS)
%       'display', false - controls plotting
%
%   See also PSTH_STATS, STIMES2BINRASTER, BINRASTER2PSTH, BINRASTER2APSTH,
%   APSTH, VIEWCELL2B, PARTITION_TRIALS and FILTERTRIALS.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   07-May-2012

%   Edit log: BH 7/5/12, 8/12/12, 8/27/12, 12/12/13, 5/20/14, 11/25/20
%             MA 09/05/25

% Default arguments
prs = inputParser;
addRequired(prs,'cellids',@(s)iscell(s)|iscellstr(s)|ischar(s)|isstring(s))
addRequired(prs,'event_type',@ischar)   % event type ('stim' or 'trial')
addRequired(prs,'event',@(s)ischar(s)|...
    (iscellstr(s)&isequal(length(s),2))|...
    isa(s,'function_handle'))   % reference event
addRequired(prs,'window',@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event, in seconds
addParameter(prs,'shevent',{}) % show event
addParameter(prs,'event_filter','none',@(s)ischar(s)|iscellstr(s))   % filter events based on properties
addParameter(prs,'filterinput',[])   % some filters need additional input
addParameter(prs,'maxtrialno',5000)   % downsample events if more than 'maxtrialno'
addParameter(prs,'data_type','real',@ischar)   % data type ('real' or 'virtual')
addParameter(prs,'first_event','',@(s)isempty(s)|ischar(s))   % exclude spikes before previous event
addParameter(prs,'last_event','',@(s)isempty(s)|ischar(s))   % exclude spikes after following events
addParameter(prs,'spike_def','',@(s)isempty(s)|ischar(s))   % use bursts or single spikes
addParameter(prs,'dt',0.001,@isnumeric)   % time resolution of the binraster, in seconds
addParameter(prs,'sigma',0.02,@isnumeric)     % smoothing kernel for the smoothed PSTH
addParameter(prs,'margin',[-0.1 0.1])  % margins for PSTH calculation to get rid of edge effect due to smoothing
addParameter(prs,'parts','all')   % partition trials
addParameter(prs,'isadaptive',false,@(s)islogical(s)|ismember(s,[0 1 2]))   % use adaptive PSTH algorithm
addParameter(prs,'forcesmoothedstat',false,@islogical)   % use smoothed PSTH for PSTH STATS, even if adaptive
addParameter(prs,'baselinewin',[-1.6 0],@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event for stat. testing, in seconds
addParameter(prs,'testwin',[0.2 1],@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event for stat. testing, in seconds
addParameter(prs,'relative_threshold',0.01,@(s)isnumeric(s)&s>=-1&s<=1)   % threshold used to assess interval limits in PSTH_STATS; negative thresholds selects the full window
addParameter(prs,'resdir',[])
addParameter(prs,'issave',true)
addParameter(prs,'display',true,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying rasters and PSTHs
addParameter(prs,'dostats',true)
parse(prs,cellids,event_type,event,window,varargin{:})
g = prs.Results;
if nargout > 2   % statistics will only be calculted if asked for
    g.dostats = true;
else
    g.dostats = false;
end

if ischar(cellids)
    cellids = {cellids};  % one cell ID
end
time = (window(1)+g.margin(1)):g.dt:(window(2)+g.margin(2));

% Directories
global DATAPATH
if isempty(g.resdir)
    g.resdir = fullfile(DATAPATH,'PSTH_results');  % results directory
    if ~isfolder(g.resdir)
        mkdir(g.resdir)
    end
end
fnmm = 'rasterpsth_.mat';   % filename for saving the result matrices

% Cell loop 
wb = waitbar(0,'Please wait...','Name','Running ultimatePSTH...');  % progress indicator
global WB
WB(end+1) = wb;
numCells = length(cellids);
numTimeBins = numel(time);
if strcmp(g.parts, '#Outcome')
    spsth_all={};
else 
    spsth_all = deal(nan(numCells,numTimeBins-200));
end

for iC = 1:numCells   % loop through the cells
    cellid = cellids{iC};
    % Load event structure
    event_type = lower(event_type(1:4));
            % Load trial events
    try
        try
            VE = loadcb(cellid,'TrialEvents');   % load events
            if strcmp(g.spike_def, 'burst')
                VS = loadcb(cellid,'BURSTSPIKES');   % load prealigned spikes
            elseif  strcmp(g.spike_def, 'single')
                VS = loadcb(cellid,'SINGLESPIKES');   % load prealigned spikes
            else
                if strcmp(g.data_type,'real')
                    VS = loadcb(cellid,'EVENTSPIKES');   % load prealigned spikes
                elseif strcmp(g.data_type,'virtual')
                    VS = loadcb(cellid,'EVENTVSPIKES');   % load prealigned virtual spikes
                else
                    error('MATLAB:CellBase:ultimate_psth','Unsupported data type.')
                end
            end
        catch %ME
            %disp('There was no behavioral protocol for this session.')
            %error(ME.message)
            fprintf('Skipping cell %s: Files not found.\n', cellid);
            continue;
        end
    catch
        error('Input argument ''event_type'' should be ''trial''.')
    end
    
    % Events, time, valid trials
    if iscell(event)   % different ref. event for test and baseline window
        event1 = event{1};   % baseline event
        event2 = event{2};   % test event
 
        event_pos1 = findcellstr(VS.events(:,1),event1);
        event_pos2 = findcellstr(VS.events(:,1),event2);
        if event_pos1 * event_pos2 == 0
            error('Event name not found');
        end
        stimes1 = VS.event_stimes{event_pos1};   % baseline spike times
        stimes2 = VS.event_stimes{event_pos2};   % test spike times
        triggerevent1 = VS.events{event_pos1,2};   % trigger event for baseline (event name may differ)
        triggerevent2 = VS.events{event_pos2,2};   % trigger event for test (event name may differ)
        valid_trials1 = filterTrials_wm(cellid,'event_type',event_type,'event',event1,...
            'event_filter',g.event_filter,'filterinput',g.filterinput);   % filter trials
        valid_trials2 = filterTrials_wm(cellid,'event_type',event_type,'event',event2,...
            'event_filter',g.event_filter,'filterinput',g.filterinput);   % filter trials
        if ~isequal(valid_trials1,valid_trials2)
            error('Valid trials should be the same for baseline and test period.')
        else
            valid_trials = valid_trials1;
        end
        [stimes1 starttimes1 endtimes1] = dynamicSpikeWindow(stimes1,VE,triggerevent1,g.first_event,g.last_event);   % restrict between previous and following event
        [stimes2 starttimes2 endtimes2] = dynamicSpikeWindow(stimes2,VE,triggerevent2,g.first_event,g.last_event);   % restrict between previous and following event
        
    else
        if isa(event,'function_handle')
            event = feval(event,cellid);   % dynamic event definition
        end
        try event_pos = findcellstr(VS.events(:,1),event);
        catch 
            fprintf('could not process from this cell')
            continue;
        end
        if event_pos == 0
            error('Event name not found');
        end
        stimes = VS.event_stimes{event_pos};
        triggerevent = VS.events{event_pos,2};   % trigger event (event name may differ)
        valid_trials = filterTrials_wm(cellid,'event_type',event_type,'event',event,...
            'event_filter',g.event_filter,'filterinput',g.filterinput);   % filter trials
             
        [stimes starttimes endtimes] = dynamicSpikeWindow(stimes,VE,triggerevent,g.first_event,g.last_event);   % restrict between previous and following event
    end

    if ~isempty(g.shevent)&&g.display==false
        EventTimes = trialevents2relativetime(VE,g.event,g.shevent);
        if strcmp(g.parts,'#Outcome')
            EventTimesR{iC}=EventTimes(1, valid_trials);
            EventTimesP{iC}=EventTimes(2,valid_trials);
            EventTimesEnd{iC}=EventTimes(3,valid_trials);
        else 
            EventTimes=EventTimes(valid_trials);
        end
    else
        EventTimes=[];
    end

    if isempty(valid_trials)
        fprintf('Skipping cell %d: No valid trials in this session.\n', cellid);
        continue;
    end

    % Calculate bin rasters
    if iscell(event)   % different ref. event for test and baseline window
        spt1 = stimes2binraster(stimes1,time,g.dt);
        spt1 = nanpadspt(time,spt1,starttimes1,endtimes1);  % replace zeros with NaNs outside the dynamic window
        spt2 = stimes2binraster(stimes2,time,g.dt);
        spt2 = nanpadspt(time,spt2,starttimes2,endtimes2);  % replace zeros with NaNs outside the dynamic window
        spt = [spt1(:,time<0) spt2(:,time>=0)];   % merge baseline and test raster
    else
        spt = stimes2binraster(stimes,time,g.dt);
        spt = nanpadspt(time,spt,starttimes,endtimes);  % replace zeros with NaNs outside the dynamic window
    end
    
    % Partition trials
    [COMPTRIALS, tags] = partition_trials(VE,g.parts);
    
    % PSTH
    switch g.isadaptive
        case {0,false}
            [psth, spsth, spsth_se] = binraster2psth(spt,g.dt,g.sigma,COMPTRIALS,valid_trials);
        case {1, true}
            [psth, spsth, spsth_se] = binraster2apsth(spt,g.dt,g.sigma,COMPTRIALS,valid_trials);
        case 2
            [psth, spsth, spsth_se] = binraster2dapsth(spt,g.dt,g.sigma,COMPTRIALS,valid_trials);
    end
    stm0 = abs(window(1)+g.margin(1)) * (1 / g.dt);   % zero point (diveding directly with 'g.dt' can result in numeric deviation from the desired integer result)
    stm = round(stm0);   % still numeric issues
    if abs(stm-stm0) > 1e-10
        error('Zero point is not an integer.')
    end
    inx = (stm+1+window(1)/g.dt):(stm+1+window(2)/g.dt);     % indices for cutting margin
    psth = psth(:,inx);
    spsth = spsth(:,inx);
    spsth_se = spsth_se(:,inx);
    
    NumPartitions = size(psth,1);
    if NumPartitions > 1   % partitions
        pspt = spt;
        spt = cell(1,NumPartitions);
        for iP = 1:NumPartitions
            spt{iP} = pspt(intersect(valid_trials,COMPTRIALS{iP}),inx);
            spsth_all{iC,iP} = spsth(iP, :);
        end
    else
        spt = spt(valid_trials,inx);
        spsth_all(iC,:) = spsth;
    end
    

    % Output statistics
    if g.dostats
        switch g.isadaptive
            case {0,false}
                statpsth = spsth;   % use smoothed PSTH for finding activation and inhibition windows in psth_stats
            case {1,2,true}
                if ~g.forcesmoothedstat
                    statpsth = psth;   % use adaptive Spike Density Function for finding activation and inhibition windows in psth_stats
                else
                    statpsth = spsth;   % use smoothed PSTH even in the adaptive cases if forced
                end
        end
        if NumPartitions > 1
            stats = cell(1,NumPartitions);   % return multiple stats if partitioning
            for iP = 1:NumPartitions
                sts = psth_stats_wm(spt{iP} ,statpsth(iP,:),g.dt,window,...
                    'baselinewin',g.baselinewin,'testwin',g.testwin,'display',g.display,...
                    'relative_threshold',g.relative_threshold);
                stats{iP} = sts;
                stats{iP}.tag = tags{iP};   % partition tag
                if g.display
                    stats{iP}.figure_handle = gcf;   % store figure handles
                    stats{iP}.axis_handle = gca;
                end
                st=stats;
            end
        else
            stats = psth_stats_wm(spt,statpsth,g.dt,window,...
                'baselinewin',g.baselinewin,'testwin',g.testwin,'display',g.display,...
                'relative_threshold',g.relative_threshold);
            if g.display
                stats.figure_handle = gcf;
                stats.axis_handle = gca;
            end
            st=stats;
        end
        if g.dostats
            stats_all(iC,:)=st;
        end
    end
    
    if g.issave && g.display   % save figure
        ncl = regexprep(cellid,'\.','_');
        fnm = ['PSTH_' ncl '.fig'];
        saveas(gcf,fullfile(g.resdir,fnm))   
        fnm2 = ['PSTH_' ncl '.jpg'];
        saveas(gcf,fullfile(g.resdir,fnm2))   
    end

    if not(numCells==1) || g.display == 0
        disp(['Cell #' num2str(iC) ' / ' num2str(numCells) ' done......'])
        close(gcf);
    end
    waitbar(iC/numCells)
end
    % Save
    
    if g.issave
        if isequal(mod(iC,50),0)   % save after every 50 pairs to prevent data loss
            save(fullfile(g.resdir,fnm),'cellids','spsth_all','time')
            disp('Autosave done.')
        end
    end
    time=time(:,inx);
    if ~isempty(g.shevent)&& g.display==false
        if strcmp(g.parts,'#Outcome')
            EventTimes = struct('Reward',EventTimesR, 'Punishment', EventTimesP, 'TrialEnds', EventTimesEnd);
        end
    end
    % Save
     if g.issave
         save(fullfile(g.resdir,fnmm),'cellids','spsth_all','time','EventTimes')
     end
    
close(wb)   % eliminate progress indicator
  
end

% -------------------------------------------------------------------------
function spt = nanpadspt(time,spt,starttimes,endtimes)

% For variable windows, change padding to NaN to ensure correct averaging
NUMtrials = size(spt,1);   % number of trials
for iT = 1:NUMtrials    % loop through trials
    inx = time < starttimes(iT);
    spt(iT,inx) = NaN;   % NaN trials before previous event
end
for iT = 1:NUMtrials    % loop through trials
    inx = time > endtimes(iT);
    spt(iT,inx) = NaN;   % NaN trials after following event
end
end