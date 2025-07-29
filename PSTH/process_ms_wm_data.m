function cleaned_data = process_ms_wm_data(resdir, varargin)
%PROCESS_MS_WM_DATA - Process and clean electrophysiology data.
%     PROCESS_MS_WM_DATA(RESDIR) performs full pipeline processing including:
%   - Spike sorting and filtering
%   - Peri-Stimulus Time Histogram (PSTH) computation
%   - Delay response categorization
%   - Theta index calculation
%   - Saving intermediate results for resuming interrupted runs
% The results are saved in a struct and includes: Normalized PSTHs, Delay 
% response matrices, Theta indices, Stats and cell IDs and Time vector.
%
% Inputs:
%   resdir - Root directory where results will be saved
%   varargin - Optional name-value pairs:
%       'spsthCtrl' - Spike data for control group (optional)
%       'spsthExp' - Spike data for experimental group (optional)
%       'stats_Ctrl' - Statistics for control group (optional)
%       'stats_Exp' - Statistics for experimental group (optional)
%       'ThetaIndexCtrl' - Theta indices for control group (optional)
%       'ThetaIndexExp' - Theta indices for experimental group (optional)
%       'window' - Time window for analysis (default: [-2, 6] seconds relative to event)
%
% Output:
%   cleaned_data - A structure containing:
%       .spsth_data - Raw and normalized PSTH data
%       .delay_response - Categorized delay activity
%       .ThetaIndex - Theta oscillation metrics
%       .stats - Statistics for delay period
%       .cellids - Cell ids
%       .time - Time vector used in PSTHs
%
% Example:
%   % Basic call with default settings
%   cleaned_data = process_ms_wm_data('E:\WorkingMemory\Data\Processed');
%
%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025

    % Parse optional inputs
    opt_args = inputParser;
    addOptional(opt_args, 'spsthCtrl', []);
    addOptional(opt_args, 'spsthExp', []);
    addOptional(opt_args, 'stats_Ctrl', []);
    addOptional(opt_args, 'stats_Exp', []);
    addOptional(opt_args, 'ThetaIndexCtrl', []);
    addOptional(opt_args, 'ThetaIndexExp', []);
    addOptional(opt_args, 'window', [-2 6]);
    parse(opt_args, varargin{:});
    
    % Unpack parsed inputs
    spsthCtrl = opt_args.Results.spsthCtrl;
    spsthExp = opt_args.Results.spsthExp;
    stats_Ctrl = opt_args.Results.stats_Ctrl;
    stats_Exp = opt_args.Results.stats_Exp;
    ThetaIndexCtrl = opt_args.Results.ThetaIndexCtrl;
    ThetaIndexExp = opt_args.Results.ThetaIndexExp;
    window = opt_args.Results.window;
    
    % Create directories
    dirs = create_result_directories(resdir);
    
    % Start Timer
    tic
    
    % Process Experimental Group : SPSTH, delay response sorting, ACG, &
    % theta index calculation
    [exp_data, cellids_exp_cl, cellids_exp] = process_group('MS_WM_EXP_cellbase', window, dirs.PSTH_Experimental, ...
        spsthExp, stats_Exp, ThetaIndexExp, resdir);
    
    % Process Control Group : SPSTH, delay response sorting, ACG, &
    % theta index calculation
    [ctrl_data, cellids_ctrl_cl, cellids_ctrl] = process_group('MS_WM_CTRL_cellbase', window, dirs.PSTH_Control, ...
        spsthCtrl, stats_Ctrl, ThetaIndexCtrl, resdir);
    
    % Combine Results
    cleaned_data = struct(...
        'spsth_data', struct(...
            'normPSTH_ctrl_cl', {ctrl_data.normPSTH_cl}, ...
            'normPSTH_exp_cl', {exp_data.normPSTH_cl}, ...
            'normPSTH_ctrl', {ctrl_data.normPSTH}, ...
            'normPSTH_exp', {exp_data.normPSTH}, ...
            'spsthCtrl_raw', {ctrl_data.spsth}, ...
            'spsthExp_raw', {exp_data.spsth} ...
        ), ...
        'delay_response', struct(...
            'ResponseCategoriExp', {exp_data.delay_response_cl}, ...
            'ResponseCategoriCtrl', {ctrl_data.delay_response_cl}, ...
            'CtrlMatrix_raw', {ctrl_data.delay_response}, ...
            'ExpMatrix_raw', {exp_data.delay_response} ...
        ), ...
        'ThetaIndex', struct(...
            'ThetaIndexCtrl', {ctrl_data.theta_index}, ...
            'ThetaIndexExp', {exp_data.theta_index} ...
        ), ...
        'stats', struct(...
            'stats_Ctrl_cl', {ctrl_data.stats_cl}, ...
            'stats_Exp_cl', {exp_data.stats_cl}, ...
            'stats_Ctrl', {ctrl_data.stats}, ...
            'stats_Exp', {exp_data.stats} ...
        ), ...
        'cellids', struct(...
            'cellids_ctrl_cl', {cellids_ctrl_cl}, ...
            'cellids_exp_cl', {cellids_exp_cl}, ...
            'cellids_ctrl', {cellids_ctrl}, ...
            'cellids_exp', {cellids_exp} ...
        ), ...
        'time', exp_data.time ...
    );
    
    % Final save
    save(fullfile(resdir, 'MS_WM_data_filtered.mat'), 'cleaned_data', '-mat');
    toc
end

function dirs = create_result_directories(resdir)
% Creates result directories structure

    dirs.ACG_Control = fullfile(resdir, 'ACG_Control');
    dirs.PSTH_Control = fullfile(resdir, 'PSTH_Control');
    dirs.ACG_Experimental = fullfile(resdir, 'ACG_Experimental');
    dirs.PSTH_Experimental = fullfile(resdir, 'PSTH_Experimental');

    for dirName = fieldnames(dirs)'
        if ~exist(dirs.(dirName{1}), 'dir')
            mkdir(dirs.(dirName{1}));
        end
    end
end

function [group_data, cellids_cl, cellids_all] = process_group(cellbase_name, window, psth_dir, ...
                                                            spsth_input, stats_input, theta_input, resdir)
% Process Exp or Ctrl data : PSTH, normalization of PSTH, delay response sorting, ACG 

    if strcmp(cellbase_name, 'MS_WM_EXP_cellbase')
        group_name = 'Exp';
    else 
        group_name = 'Ctrl';
    end
    
    [success, group_data, cellids_cl, cellids_all, step_suffix] = try_resume(resdir, group_name); % Try to resume from autosaved or previously saved files
    if success
        % Skip steps already done
        switch step_suffix
            case 'theta_done'
                return;  % All processing done
            case 'response_done'
                spsth_input = group_data.spsth;
                fprintf('Resuming after response sorter...\n');
            case 'psth_done'
                spsth_input = group_data.spsth;
                fprintf('Resuming after PSTH computation...\n');
        end
    else
        % Start fresh
        group_data = struct();
    end
    
    choosecb(cellbase_name);
    loadcb; % Load TheMatrix, CELLIDLIST, ANALYSES
    cellids_all = CELLIDLIST';
    
    % Select well-isolated units in the MS
    cellids=selectcell('"ID">20&"Lratio"<0.15&ismember("Anatomy",''MS'')');
    
    % Run PSTH loop if not provided or not resumed
    if isempty(spsth_input) && ~isfield(group_data, 'spsth')
    fprintf('Processing %s: Running ultimate_psth_wm...\n', group_name);
    
        % Calculate smoothed PSTH, stats for delay period & save PSTH figures
        [spsth, stats_output, time, ~] = ultimate_psth_wm(cellids, ...
            'trial', 'FixationBeginning', window, 'resdir', psth_dir, ...
            'display', true, 'dt', 0.001, 'sigma', 0.02, 'isadaptive', 0, ...
            'maxtrialno', Inf, 'relative_threshold', 0.01, 'event_filter', 'lowfixation_wm');
        
        group_data.spsth = spsth;
        group_data.stats = stats_output;
        group_data.time = time;
        
        autosave_step(resdir, group_name, 'psth_done', locals2struct());
    
    elseif isfield(group_data, 'spsth')
        % Already loaded via resume
    else
        % Use input values
        group_data.spsth = spsth_input;
        group_data.stats = stats_input;
        group_data.time = linspace(window(1), window(2), 8001);
    end
    
    % Normalize and clean PSTH data
    [~, normPSTH, ~, ~]  = normalize_psth(group_data.spsth, group_data.time);
    if isa(normPSTH, 'double') && size(normPSTH, 1) > 1
        valid_idx = sum(isnan(normPSTH), 2) == 0;
    else
        valid_idx = true(size(normPSTH, 1), 1); % fallback
    end
    
    group_data.normPSTH = normPSTH;
    group_data.normPSTH_cl = normPSTH(valid_idx, :);
    group_data.stats_cl = group_data.stats(valid_idx, :);
    cellids_cl = cellids(valid_idx)';
    
    % Run delay response sorter if not resumed
    if ~isfield(group_data, 'delay_response') 
        if findprop('delay_response')~=0
                delay_response = getvalue('delay_response', cellids);
        else 
            fprintf('Processing %s: Running wm_responsesorter...\n', group_name);
            
            % Group cells by response during the delay period
            wm_responsesorter(CELLIDLIST, 'true', window, psth_dir, cellbase_name);
            delay_response = getvalue('delay_response', cellids);
        end
        group_data.delay_response_cl = delay_response(valid_idx);
        group_data.delay_response = delay_response;
        
        autosave_step(resdir, group_name, 'response_done', locals2struct());
    end
    
    % Compute Theta Index if not resumed
    if ~isfield(group_data, 'theta_index') &&   isempty(theta_input)
    
        fprintf('Processing %s: Running acgmod...\n', group_name);

        % Produce ACGs, save figs & calculate theta index
        [~, ~, theta_index] = acgmod(cellids_cl, ...
            1, 'resdir', strrep(psth_dir, 'PSTH', 'ACG'), 'dt', 0.001);
    
        group_data.theta_index = theta_index;
    
        autosave_step(resdir, group_name, 'theta_done', locals2struct());
    end
end

function [success, group_data, cellids_cl, cellids_all, step_suffix] = try_resume(resdir, group_name)
% Try to resume from autosaved file.

    % Initialize all outputs to default values
    group_data = struct();
    cellids_cl = [];
    cellids_all = [];

    stages = {'psth_done', 'response_done', 'theta_done'};
    latest_file = '';
    latest_time = 0;

    % Look for files corresponding to this group only
    for i = 1:length(stages)
        fname = fullfile(resdir, ['MS_WM_data_filtered_', group_name, '_', stages{i}, '.mat']);
        if exist(fname, 'file')
            fileInfo = dir(fname);
            modifiedTime = datenum(fileInfo.date); % convert date string to numeric
            if modifiedTime > latest_time
                latest_time = modifiedTime;
                latest_file = fname;
            end
        end
    end

    if isempty(latest_file)
        success = false;
        step_suffix = '';
        return;
    end

    % Extract stage (step_suffix) from filename
    [~, name, ~] = fileparts(latest_file);
    parts = split(name, '_');
    step_suffix = strjoin(parts(end-1:end), '_');  

    % Load variables from file
    fprintf('Resuming %s from: %s\n', group_name, latest_file);
    try
        loaded = load(latest_file);

        % Assign outputs only if they exist in the file
        group_data = loaded.group_data;
        if isfield(loaded, 'cellids_cl')
            cellids_cl = loaded.cellids_cl;
        else
            cellids_cl = [];
        end
        
        if isfield(loaded, 'cellids_all')
            cellids_all = loaded.cellids_all;
        else
            cellids_all = [];
        end

        success = true;
        fprintf('Successfully resumed from stage: %s\n', step_suffix);

    catch ME
        warning('Failed to load resume file: %s', latest_file);
        disp(ME.message);
        success = false;
        step_suffix = '';
    end
end

function s = locals2struct()
    % List of variables to save - include all expected ones
    varnames = {'group_data', 'cellids_cl', 'cellids_all'};
    
    % Build structure from caller workspace
    s = struct();
    for i = 1:length(varnames)
        varname = varnames{i};
        if evalin('caller', ['exist(''' varname ''', ''var'')']) == 1
            s.(varname) = evalin('caller', varname);
        else
            warning('Variable "%s" does not exist in caller workspace.', varname);
        end
    end
end

function autosave_step(resdir, group_name, step_suffix, data)
% Autosave to not lose progress
    filename = fullfile(resdir, ['MS_WM_data_filtered_', group_name, '_', step_suffix, '.mat']);
    save(filename, '-struct', 'data');
    fprintf('Saved intermediate result: %s\n', filename);
end