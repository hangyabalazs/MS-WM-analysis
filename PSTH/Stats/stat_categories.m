function [meansGroups, SEGroups, WpGroups, meansGroupsSub, SEGroupsSub, WpGroupsSub] = stat_categories(cleaned_data, resdir)
%STAT_CATEGORIES Statistical summaries and plotting of firing rates.
%   [MEANSGROUPS, SEGROUPS, WPGROUPS, MEANSGROUPSSUB, SEGROUPSSUB, WPGROUPSSUB] = ...
%       STAT_CATEGORIES(CLEANED_DATA, RESDIR) computes statistical summaries
%   and generates plots comparing firing rates across task phases and response
%   categories in experimental and control groups.
%
%   Inputs:
%       cleaned_data : Struct containing preprocessed PSTH data and response
%                      categories.
%       resdir       : Directory path where figures are saved and feedback data
%                      is loaded from.
%
%   Outputs:
%       meansGroups    : Mean firing rates per group (normal data).
%       SEGroups       : Standard errors per group (normal data).
%       WpGroups       : Statistical test results; p values (normal data).
%       meansGroupsSub : Mean firing rates after baseline subtraction.
%       SEGroupsSub    : Standard errors after baseline subtraction.
%       WpGroupsSub    : Statistical test results after baseline subtraction.
%
%   Example:
%       [means, SE, stats, meansSub, SESub, statsSub] = stat_categories(dataStruct, 'results/');
%
%   See also MEANBARGROUPS, SPSTH_PARTS.
%
%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025

    close all; % Avoid confusion
    
    % Define time window & resolution
    g.window = [-2 6]; % Time window for analysis in seconds
    g.dt = 0.001;      % Time resolution (sampling interval)

    % Significance level
    alpha = 0.05;
    
    % Extract raw PSTHs and normalized PSTHs
    spsthExp = cleaned_data.spsth_data.spsthExp_raw;
    spsthCtrl = cleaned_data.spsth_data.spsthCtrl_raw;
    normPSTH_exp = cleaned_data.spsth_data.normPSTH_exp;
    normPSTH_ctrl = cleaned_data.spsth_data.normPSTH_ctrl;
    
    % Remove entries with NaNs in normalized PSTHs
    spsthExp = spsthExp(sum(isnan(normPSTH_exp), 2) == 0, :);
    spsthCtrl = spsthCtrl(sum(isnan(normPSTH_ctrl), 2) == 0, :);
    
    % Response categories for experimental and control groups
    ME = categorical(cleaned_data.delay_response.ResponseCategoriExp);
    MC = categorical(cleaned_data.delay_response.ResponseCategoriCtrl);
    
    % Define group names
    group_names = {'Inh', 'Act', 'Inh-Act', 'Act-Inh', 'NonResp'};
    
    % Create time vector and delay indices for analysis windows
    time_vector = g.window(1):g.dt:g.window(2);
    delay_start = 0.2;
    delay_end = 1;
    delay_inx = find(time_vector >= delay_start & time_vector <= delay_end);
    
    % Prepare paths and load feedback firing rate data (reward and punishment)
    datapath = fullfile(resdir, 'feedback');
    if ~exist(datapath, 'dir')
        error('Folder %s does not exist.', datapath);
    end
    
    filename = fullfile(datapath, 'meanFR.mat');
    
    if exist(filename, 'file')
        disp(['Loading existing data: ' filename]);
        load(filename, 'meanFR');
        feedback_meanFR = meanFR;
    else
        disp(['No data found. Running spsth_parts() and saving to: ' filename]);
        feedback_meanFR = spsth_parts(datapath);
    end
    
    % Process reward and punishment data for Ctrl and Exp groups
    [mn_rewCtrl, mn_punishCtrl] = process_feedback_data(feedback_meanFR, normPSTH_ctrl, 'Ctrl');
    [mn_rewExp, mn_punishExp] = process_feedback_data(feedback_meanFR, normPSTH_exp, 'Exp');
    
    % Define firing rate conditions for statistical comparison
    conditions = {
        'Overall', mean(spsthExp, 2), mean(spsthCtrl, 2);
        'Baseline', mean(spsthExp(:, time_vector >= -1.6 & time_vector <= 0), 2), mean(spsthCtrl(:, time_vector >= -1.6 & time_vector <= 0), 2);
        'Tone', mean(spsthExp(:, time_vector >= 0 & time_vector <= delay_start), 2), mean(spsthCtrl(:, time_vector >= 0 & time_vector <= delay_start), 2);
        'Delay', mean(spsthExp(:, delay_inx), 2), mean(spsthCtrl(:, delay_inx), 2);
        'Reward', mn_rewExp, mn_rewCtrl;
        'Punishment', mn_punishExp, mn_punishCtrl;
    };
    
    group_data = cell(6, 10); % Preallocate container for grouped data
    
    % Organize data by category and condition for both groups
    for iC = 1:6
        mnExp = conditions{iC, 2};
        mnCtrl = conditions{iC, 3};
        for group = 1:length(group_names)
            groupIndicesE = (ME == group_names{group});
            groupIndicesC = (MC == group_names{group});
            for iP = 1:10
                % Index column alternates Ctrl (odd) and Exp (even)
                columnIndex = (group - 1) * 2 + mod(iP, 2) + 1;
                if mod(iP, 2) == 0
                    group_data{iC, columnIndex} = mnExp(groupIndicesE, :);
                else
                    group_data{iC, columnIndex} = mnCtrl(groupIndicesC, :);
                end
            end
        end
    end
    
    % Define group names
    groupAll = {'Inh Exp','Inh Ctrl','Act Exp','Act Ctrl','Inh-Act Exp','Inh-Act Ctrl','Act-Inh Exp','Act-Inh Ctrl','NonResp Exp','NonResp Ctrl'};
    
    meansGroups = cell(6,1);
    SEGroups = cell(6,1);
    WpGroups = cell(6,1);
    
    % Generate bar graph, calculate group means, standard errors, and perform stats for
    % non-transformed data
    for iC = 1:6
        cond = group_data(iC, :);
        [axes_handles{iC}, ~, WpGroups{iC}, meansGroups{iC}, SEGroups{iC}] = meanbargroups(cond, groupAll, alpha, conditions{iC}, 0); % Save bar graph axes and stats
    end
    
    % Baseline subtraction for Exp and Ctrl groups
    baselineExp = mean(spsthExp(:, time_vector >= -1.6 & time_vector <= 0), 2);
    baselineCtrl = mean(spsthCtrl(:, time_vector >= -1.6 & time_vector <= 0), 2);
    
    conditionsSub = {
        'Tone', mean(spsthExp(:, time_vector >= 0 & time_vector <= delay_start), 2) - baselineExp, mean(spsthCtrl(:, time_vector >= 0 & time_vector <= delay_start), 2) - baselineCtrl;
        'Delay', mean(spsthExp(:, delay_inx), 2) - baselineExp, mean(spsthCtrl(:, delay_inx), 2) - baselineCtrl;
        'Reward', mn_rewExp - baselineExp, mn_rewCtrl - baselineCtrl;
        'Punishment', mn_punishExp - baselineExp, mn_punishCtrl - baselineCtrl;
    };
    
    group_dataSub = cell(4, 10);
    
    % Organize baseline-subtracted data by category and condition
    for iC = 1:4
        mnExp = conditionsSub{iC, 2};
        mnCtrl = conditionsSub{iC, 3};
        for group = 1:length(group_names)
            groupIndicesE = (ME == group_names{group});
            groupIndicesC = (MC == group_names{group});
            for iP = 1:10
                columnIndex = (group - 1) * 2 + mod(iP, 2) + 1;
                if mod(iP, 2) == 0
                    group_dataSub{iC, columnIndex} = mnExp(groupIndicesE, :);
                else
                    group_dataSub{iC, columnIndex} = mnCtrl(groupIndicesC, :);
                end
            end
        end
    end
    
    meansGroupsSub = cell(4,1);
    SEGroupsSub = cell(4,1);
    WpGroupsSub = cell(4,1);
    
    % Generate bar graph and calculate means, SEs, and stats for baseline-subtracted data
    for iC = 1:4
        cond = group_dataSub(iC, :);
        [axes_handles{iC+6}, ~ , WpGroupsSub{iC}, meansGroupsSub{iC}, SEGroupsSub{iC}] = meanbargroups(cond, groupAll, alpha, conditionsSub{iC}, 0); % Save bar graph axes and stats
    end
    
    % Load psth partitioned by outcome or run analysis again if unavailable
    datapath = fullfile(resdir,'partitionedbyoutcome');
    try
        % Load Experimental Data
        exp_file = fullfile(datapath,'Exp\rasterpsth_.mat');
        load(exp_file, 'spsth_all');
        spsthExpRew = cell2mat(spsth_all(:,1));  
        spsthExpPun = cell2mat(spsth_all(:,2));
        
        % Load Control Data
        ctrl_file = fullfile(datapath,'Ctrl\rasterpsth_.mat');
        load(ctrl_file, 'spsth_all');
        spsthCtrlRew = cell2mat(spsth_all(:,1));
        spsthCtrlPun = cell2mat(spsth_all(:,2));
        
        clear spsth_all;
    catch
        % Ask the user if they want to rerun the analysis
        userResponse = input('No results found for partitioned by outcome psth analysis. Run analysis? (y/n): ', 's');
        if strcmpi(userResponse, 'y')
            run_psth_byoutcome(datapath);
        else 
            error('No analysis files found. Cannot proceed.')
        end
        
        % Load Experimental Data
        exp_file = fullfile(datapath,'Exp\rasterpsth_.mat');
        load(exp_file, 'spsth_all');
        spsthExpRew = cell2mat(spsth_all(:,1));  
        spsthExpPun = cell2mat(spsth_all(:,2));
        
        % Load Control Data
        ctrl_file = fullfile(datapath,'Ctrl\rasterpsth_.mat');
        load(ctrl_file, 'spsth_all');
        spsthCtrlRew = cell2mat(spsth_all(:,1));
        spsthCtrlPun = cell2mat(spsth_all(:,2));
        
        clear spsth_all;
    end
    
    % Subtract baseline mean from psth
    ExpRewGroups = subtract_meanbaseline(spsthExpRew, 5,  ME, group_names, time_vector); % WM/Rewarded trials
    ExpPunGroups = subtract_meanbaseline(spsthExpPun, 5,  ME, group_names, time_vector); % WM/Punished trials
    CtrlRewGroups = subtract_meanbaseline(spsthCtrlRew, 5, MC, group_names, time_vector); % Ctrl/Rewarded trials
    CtrlPunGroups = subtract_meanbaseline(spsthCtrlPun, 5, MC, group_names, time_vector); % Ctrl/Punished trials
    
    % Organize psth data for Exp & Ctrl into one cell array 
    Exp_allgroups = organize_data(5, ExpRewGroups, ExpPunGroups, time_vector);
    Ctrl_allgroups = organize_data(5, CtrlRewGroups, CtrlPunGroups, time_vector);
    
    % Compare Correct vs Incorrect
    [axes_handles{11}, ~, Wp1, means1, SEs1] = meanbargroups(Exp_allgroups, groupAll, alpha, 'DR-2AFC WM - Delay', 1); % WM 
    [axes_handles{12}, ~, Wp2, means2, SEs2] = meanbargroups(Ctrl_allgroups, groupAll, alpha, 'DR-2AFC Control - Delay', 1); % Control
    
    % Plotting parameters for figures
    fig1_width = 15; % cm
    fig1_height = 20; % cm
    figure('Units', 'centimeters', 'Position', [0, 0, fig1_width, fig1_height]);
    t1 = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    tick_length_cm = 0.2;
    tick_length_y = tick_length_cm / fig1_width;
    tick_length_x = tick_length_cm / fig1_height;
    
    % Copy plots from meanbargroups output to tiled layout for non-transformed data
    for i = 1:12
        ax = nexttile(t1, i); % New tile
        children = get(axes_handles{i}, 'Children'); % Copy all children (lines, patches, etc) from original axes
        copyobj(children,  ax);
        ylabel(ax, get(get(axes_handles{i}, 'YLabel'), 'String')); % Copy YLabel
        title(ax, get(get(axes_handles{i}, 'Title'), 'String')); % Copy Title
        xticks(ax, get(axes_handles{i}, 'XTick')); % Copy Xticks
        xticklabels(ax, get(axes_handles{i}, 'XTickLabel')); % Copy Xticks labels
        labels = get(get(axes_handles{i}, 'Legend'),'String');
        set(ax, 'TickDir', 'out', 'Box', 'off', 'TickLength', [tick_length_x, tick_length_y]); % Set plot 
        if i >= 11
            p1 = patch(NaN, NaN, 'w', 'EdgeColor', 'green', 'LineWidth', 1);
            p2 = patch(NaN, NaN, 'w', 'EdgeColor', 'red', 'LineWidth', 1);
        else 
            p1 = patch(NaN, NaN, 'w', 'EdgeColor', 'red', 'LineWidth', 1);
            p2 = patch(NaN, NaN, 'w', 'EdgeColor', 'blue', 'LineWidth', 1);
        end
        legend([p1, p2], labels, 'Location', 'best','Box','off');
         % Create legend
        lgd = legend(ax, [p1, p2], labels, 'Location', 'best', 'Box', 'off', ...
                      'FontSize', 8);
    
        % Reduce internal spacing
        lgd.ItemTokenSize = [10 10]; % Small square legend icons
    end
    
    set(gcf, 'Renderer', 'painters');
    
    % Save figure 
    saveas(gcf, fullfile(resdir, 'Stat panel S6.jpg'));
    saveas(gcf, fullfile(resdir, 'Stat panel S6.svg'));

end

function run_psth_byoutcome(datapath)
% Run psth analysis partitioned by outcome

    if ~isfolder(datapath)
        mkdir(datapath);
    end
    % Set up dataset names and alignment conditions
    datasets = {'Exp', 'Ctrl'};
        % Loop through each dataset (Experimental and Control)
    for d = 1:length(datasets)
        if strcmp(datasets{d}, 'Exp')
            cellbase_name = 'MS_WM_EXP_cellbase';
        else
            cellbase_name = 'MS_WM_CTRL_cellbase';
        end
        choosecb(cellbase_name);
        loadcb;  % Reload database
        cellids = selectcell('"ID">20&"Lratio"<0.15&ismember("Anatomy",''MS'')');

        % Loop through Reward and Punishment alignments
            resdir1 = fullfile(datapath,  datasets{d});
            
            % Check if directory and results file exist
            datafile = fullfile(resdir1, 'rasterpsth_.mat');
            if exist(datafile, 'file')
                fprintf('Results already exist for %s-%s. Skipping analysis.\n', datasets{d});
            else
                if ~isfolder(resdir1)
                    mkdir(resdir1);
                end
                
                % Ask the user if they want to rerun the analysis
                userResponse = input(['No results found for ', datasets{d}, '-', '. Run analysis? (y/n): '], 's');
                if strcmpi(userResponse, 'y')
                    % Run ultimate_psthloop
                    fprintf('Running analysis for %s-%s...\n', datasets{d});
                    ultimate_psth_wm(cellids, 'trial', 'FixationBeginning', [-2 6], ...
                        'parts', '#Outcome', 'resdir', resdir1, 'dt', 0.001, 'sigma', 0.02, ...
                        'isadaptive', 0, 'maxtrialno', Inf, 'relative_threshold', 0.01, ...
                        'shevent', {'RewardBeginning', 'PunishmentBeginning', 'TrialEnd'},...
                        'event_filter', 'lowfixation_wm', 'display', false);
                else
                    fprintf('Skipping analysis for %s-%s.\n', datasets{d});
                    continue;
                end
            end

            % Load and process results
            load(datafile, 'spsth_all', 'time', 'EventTimes');
     end
end


function [mn_rew, mn_punish] = process_feedback_data(feedback_meanFR, normPSTH, type)
% Extract reward and punishment mean firing rates
% after removing NaN rows based on normalized PSTH data

    nan_rows = all(isnan(normPSTH), 2); % Rows with all NaNs
    
    if strcmp(type, 'Exp')
        mn_rew = feedback_meanFR{1,1}(~nan_rows, :);
        mn_punish = feedback_meanFR{1,2}(~nan_rows, :);
    else
        mn_rew = feedback_meanFR{2,1}(~nan_rows, :);
        mn_punish = feedback_meanFR{2,2}(~nan_rows, :);
    end
end

function group_data = subtract_meanbaseline(spsth_data, num_groups, grouping_matrix, Labels, time_vector)
% Subtract baseline mean from spsth data

    % Define baseline Indices
    baselineIndices = find(time_vector >= -1.6 & time_vector <= 0);
    
    % Get dimensions of spsth_data
    [numPSTH, numTimePoints] = size(spsth_data);
    subtractedmean_spsth = nan(numPSTH, numTimePoints);
    
    % Loop through each PSTH to normalize using baseline period
    for iPSTH = 1:numPSTH
        mean_baseline = mean(spsth_data(iPSTH, baselineIndices)); % Mean of baseline period
        subtractedmean_spsth(iPSTH, :) = (spsth_data(iPSTH, :) - mean_baseline) ; % Subtracted Mean SPSTH 
    end
    
    group_data = cell(1, num_groups);
    
    for group = 1:num_groups
        group_data{group} = subtractedmean_spsth(grouping_matrix == Labels{group}, :);
    end

end

function spsth_allgroups = organize_data(num_groups, spsth_groups1, spsth_groups2, time_vector)
% Organize spsth data into one cell array

    % Define time window indices; full delay, 1st half and 2nd half
    delay_start = 0.2;
    delay_end = 1;
    delay_inx = find(time_vector >= delay_start & time_vector <= delay_end);
    
    % Organize categorized spsth data into the same cell array
    spsth_data = cell(1, num_groups*2);
    
    for group = 1:num_groups
        rowIndex1 = (group - 1) * 2 + 1; 
        rowIndex2 = (group - 1) * 2 + 2;
        
        spsth_data{rowIndex1} = spsth_groups1{group};
        spsth_data{rowIndex2} = spsth_groups2{group};
    end
    
    % Compute  mean FR per catgeory 
    spsth_allgroups = cell(1, num_groups*2);
    
    for i = 1:length(spsth_data)
        spsth_allgroups{i} = mean(spsth_data{i}(:, delay_inx), 2);
    end

end