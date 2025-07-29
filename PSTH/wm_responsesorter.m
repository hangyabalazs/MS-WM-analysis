function responses = wm_responsesorter(cellids, issave, wn, response_resdir, cohort_name)
% WM_RESPONSESORTER Categorizes MS cell responses using PSTH analysis.
%
%   WM_RESPONSESORTER(CELLIDS, ISSAVE, WN, RESPONSE_RESDIR, COHORT_NAME) analyzes neural activity
%   during delay for Exp and Ctrl. It calculates non-adaptive PSTHs (see ULTIMATE_PSTH),
%   performs statistical tests (Mann-Whitney U-test), and categorizes each cell into one of 5 groups
%   based on statistically significant activity during the delay (p < 0.01):
%       'Inh'     - Only inhibition
%       'Act'     - Only activation
%       'Inh-Act' - Inhibited then activated
%       'Act-Inh' - Activated then inhibited
%       'NonResp' - No significant response
%   Adds category result to CellBase under property 'delay_response'.
%
%   INPUTS:
%       CELLIDS         - Cell ID list 
%       ISSAVE          - whether to save results to disk
%       WN              - PSTH window
%       RESPONSE_RESDIR - Directory to store output files
%       COHORT_NAME     - String identifier for saving files
%
%   OUTPUT:
%       RESPONSES       - Cell array with response categories per cell
%
%   See also: ULTIMATE_PSTH, PSTH_STATS, LRATIO
%
%   Adapted from:
%   Balazs Hangya, Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hangya.balazs@koki.mta.hu
%   11-Nov-2018

    % Initial Setup

    % Ensure result directory exists
    if ~isfolder(response_resdir)
        mkdir(response_resdir);
    end
    
    % Input argument handling
    narginchk(0,6); 
    if nargin < 2
        issave = true;   % Default saving behavior
        wn = [-2 6];     % Default window if not provided
    end
    
    % Progress bar setup
    wb = waitbar(0, 'Please wait...', 'Name', 'Running ResponseSorter...');
    numCells = length(cellids);
    
    % Load CellBase data
    choosecb(cohort_name);
    load(getpref('cellbase', 'fname'), 'CELLIDLIST');
    
    % Parameters for PSTH Analysis
    
    alignevent = 'FixationBeginning';   % Event to align trials to
    partition = 'all';                  % Trial partitioning mode
    dt = 0.001;                         % Time resolution (seconds)
    sigma = 0.02;                       % Smoothing sigma for PSTH
    bwin = [-1.6 0];                    % Baseline window for stats
    twin = [0.2 1];                     % Test window for stats
    propname = 'delay_response';        % Property name in CellBase
    
    % Create property in CellBase if not present
    if ~ismember(propname, listtag('prop'))
        insertdata([CELLIDLIST', num2cell(nan(size(CELLIDLIST')))], ...
                   'type', 'property', 'name', propname);
    end
    
    % Main Loop: Process Each Cell
    
    responses = cell(numCells, 1); 
    
    for iC = 1:numCells
        cellid = cellids{iC}; % Current cell ID
        
        try
            stats1 = rasterPSTH(cellid, alignevent, partition, wn, dt, sigma, bwin, twin);
        catch
            fprintf('Skipping cell %s: Files not found.\n', cellid);
            responses{iC} = 'Error';
            continue;
        end
    
        % Determine response type based on stats
        stat = stats1{1};
    
        if stat.inhibition_start < stat.activation_start && ...
                stat.Wpi < 0.01 && stat.Wpa < 0.01
            setvalue(cellid, propname, 'Inh-Act');
            resp = 'Inh-Act';
        elseif stat.inhibition_start > stat.activation_start && ...
                stat.Wpi < 0.01 && stat.Wpa < 0.01
            setvalue(cellid, propname, 'Act-Inh');
            resp = 'Act-Inh';
        elseif stat.Wpi < 0.01 && (stat.Wpa > 0.01 || isnan(stat.Wpa))
            setvalue(cellid, propname, 'Inh');
            resp = 'Inh';
        elseif stat.Wpa < 0.01 && (stat.Wpi > 0.01 || isnan(stat.Wpi))
            setvalue(cellid, propname, 'Act');
            resp = 'Act';
        else
            setvalue(cellid, propname, 'NonResp');
            resp = 'NonResp';
        end
    
        responses{iC} = resp;
    
        % Update progress
        disp(['Cell #' num2str(iC) ' / ' num2str(numCells) ' done......']);
        waitbar(iC / numCells, wb);
        close all; % Close any open figures
    end
    
    % Save Results
    
    close(wb); % Close progress bar
    
    % Save output file if requested
    if issave
        filename = fullfile(response_resdir, ['response_delay_', cohort_name, '.mat']);
        save(filename, 'responses');
    end
    
end 
    
function stats1 = rasterPSTH(cellid, alignevent, partition, wn, dt, sigma, bwin, twin)
% PSTH calculation wrapper
    
    TE = loadcb(cellid, 'TrialEvents');     % Load trial events
    SP = loadcb(cellid, 'EVENTSPIKES');     % Load spike times
    
    % Check for trial count consistency
    fld = fieldnames(TE);
    if ~isequal(length(SP.event_stimes{1}), length(TE.(fld{1})))
        error('MATLAB:vpresponsesorter:rasterPSTH:trialMismatch',...
            'Trial number mismatch between TrialEvents and EVENTSPIKES.')
    end
    
    % Calculate PSTH
    [ ~, stats1, ~, ~] = ultimate_psth_wm(cellid, 'trial', alignevent, wn, ...
        'dt', dt, 'sigma', sigma, 'parts', partition, 'isadaptive', 0, ...
        'maxtrialno', Inf, 'relative_threshold', 0.01, 'display', false, ...
        'event_filter', 'lowfixation_wm', 'baselinewin', bwin, 'testwin', twin,'event_filter', 'lowfixation_wm');
    
    % Wrap in cell if only one result
    if ~iscell(stats1)
        stats1 = {stats1};
    end

end 