function stat_panel_fig2(resdir, cleaned_data, varargin)
%STAT_PANEL_FIG2 Generates statistical panel figures in Fig2.
%   STAT_PANEL_FIG2(RESDIR, CLEANED_DATA) compares experimental and control data 
%   across various conditions (e.g., overall, baseline, cue, delay, reward, 
%   punishment). It computes the mean and standard error for each condition, 
%   generates cumulative distribution function (CDF) plots, and stores the results 
%   in a table saved as a CSV file. The resulting figure is saved as an SVG in the 
%   specified directory (resdir).
%
%   Inputs:
%       RESDIR         - Path to the directory where the figure and CSV file 
%                        will be saved.
%       CLEANED_DATA   - Contains the cleaned data with fields as expected 
%                        from process_ms_wm_data.m.
%
%   Optional parameter-value pairs:
%       'window'       - Time window relative to the reference event. Default: [-2 6].
%       'dt'           - Resolution of PSTH data. Default: 0.001
%       'display'      - If true, displays the plots. Default: true.
%       'save_fig'     - If true, saves the generated figure. Default: true.
%       'save_csv'     - If true, saves the results table as a CSV file. Default: true.
%
%   Example:
%       stat_panel_fig2('results/', cleaned_data, 'window', [-1 6], 'display', false);
%
% See also SPSTH_PARTS.
%
%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025

    % Default Arguments and Input Parsing
    prs = inputParser;
    addRequired(prs, 'resdir', @(x) ischar(x) || isstring(x)); % Path to save results
    addRequired(prs, 'cleaned_data', @isstruct); % Cleaned psth data structure
    addOptional(prs, 'window', [-2 6], @(x) isnumeric(x) && numel(x) == 2); % Time window
    addOptional(prs,'dt', 0.001, @isnumeric)   % time resolution of the binraster, in seconds
    addOptional(prs, 'display', true, @(x) islogical(x)); % Display flag for plots
    addOptional(prs, 'save_fig', true, @(x) islogical(x)); % Save figure flag
    addOptional(prs, 'save_csv', true, @(x) islogical(x)); % Save CSV flag
    parse(prs, resdir, cleaned_data, varargin{:});
    g = prs.Results;
    
    close all % Close open figures
    
    % Data preparation
    
    % Extract smoothed PSTH data from struct
    spsthExp = cleaned_data.spsth_data.spsthExp_raw; % Smoothed PSTH for WM task (with NaNs)
    spsthCtrl = cleaned_data.spsth_data.spsthCtrl_raw; % Smoothed PSTH for control task (with NaNs)
    normPSTH_ctrl = cleaned_data.spsth_data.normPSTH_ctrl; % Normalized PSTH for control (NaNs removed)
    normPSTH_exp = cleaned_data.spsth_data.normPSTH_exp; % Normalized PSTH for WM (NaNs removed)
    
    % Remove rows containing NaNs based on normalized PSTH
    spsthExp = spsthExp(sum(isnan(normPSTH_exp),2) == 0, :); % Keep only rows without NaNs for WM
    spsthCtrl = spsthCtrl(sum(isnan(normPSTH_ctrl),2) == 0, :); % Keep only rows without NaNs for control
    
    % Define labels for legend entries in plots
    WM_label = 'WM';
    Ctrl_label = 'Ctrl';
    
    % Define time vector based on input window and resolution
    time_vector = g.window(1):g.dt:g.window(2);
    
    % Define delay period window
    delay_start = 0.2; 
    delay_end = 1;
    delay_inx = find(time_vector >= delay_start & time_vector <= delay_end); % Indices for delay period
    
    % Reward & punishment event times and mean FR
    % Define path where feedback mean FR data is stored
    datapath = fullfile(resdir, 'feedback');
    
    % Check if the directory exists, otherwise throw an error
    if ~exist(datapath, 'dir')
        mkdir(datapath);
    end
    
    % Construct full path to the .mat file
    filename = fullfile(datapath, 'meanFR.mat');
    
    % Load existing feedback data if available; otherwise, compute it
    if exist(filename, 'file')
        disp(['Loading existing data from: ' filename]);
        load(filename, 'meanFR');
        feedback_meanFR = meanFR;
    else
        % If not found, compute feedback PSTH parts and save them
        disp(['Running spsth_parts() and saving to: ' filename]);
        feedback_meanFR = spsth_parts(datapath);  % Compute and return mean firing rates for reward and punishment events (Exp and Ctrl)
    end
    
    % Compute mean FR for reward and punishment events in control condition
    [mn_rewCtrl, mn_punishCtrl] = process_feedback_data(feedback_meanFR, normPSTH_ctrl, 'Ctrl');
    
    % Compute mean FR for reward and punishment events in experimental condition
    [mn_rewExp, mn_punishExp] = process_feedback_data(feedback_meanFR, normPSTH_exp, 'Exp');
    
    % Define list of experimental conditions with corresponding mean values
    conditions = {
        'Overall', mean(spsthExp(:, :), 2), mean(spsthCtrl(:, :), 2); % Mean over full window
        'Baseline', mean(spsthExp(:, time_vector >= -1.6 & time_vector <= 0), 2), ...
                    mean(spsthCtrl(:, time_vector >= -1.6 & time_vector <= 0), 2); % Pre-cue baseline
        'Tone', mean(spsthExp(:, time_vector >= 0 & time_vector <= delay_start), 2), ...
                mean(spsthCtrl(:, time_vector >= 0 & time_vector <= delay_start), 2); % Cue period
        'Delay', mean(spsthExp(:, delay_inx), 2), ...
                 mean(spsthCtrl(:, delay_inx), 2); % Delay period
        'Reward', mn_rewExp, mn_rewCtrl; % Feedback-related: reward
        'Punishment', mn_punishExp, mn_punishCtrl; % Feedback-related: punishment
    };
    
    % Get valid cell IDs (already aligned with spsthExp and spsthCtrl)
    cellids_exp_valid = cleaned_data.cellids.cellids_exp_cl;
    cellids_ctrl_valid = cleaned_data.cellids.cellids_ctrl_cl;
    
    % Initialize table to store per-mouse baseline FR
    MouseTable = table([], [], [], [], [], ...
        'VariableNames',{'MouseID','Group','Condition','MeanFR','SE_FR'},'RowNames', {});
    
    % Mouse identifiers
    mice_exp_list = {'NWM15', 'NWM5', 'NWM9', 'NWM7'};
    mice_ctrl_list = {'NWM21', 'NWM19', 'NWM22'};
    
    % Assign unique colors
    colors_exp = parula(length(mice_exp_list)); % One color per Exp mouse
    colors_ctrl = cool(length(mice_ctrl_list)); % One color per Ctrl mouse
    
    for condIdx = 1:size(conditions, 1)
        % Get current condition name and data
        condName = conditions{condIdx, 1};
        expData = conditions{condIdx, 2};  % Neuron x 1 vector of mean FRs for Exp
        ctrlData = conditions{condIdx, 3};  % Neuron x 1 vector of mean FRs for Ctrl
        
        % Process Experimental Mice
        for i = 1:length(mice_exp_list)
            mouse_id = mice_exp_list{i};
            
            % Find indices where cellid starts with mouse_id
            idx = cellfun(@(x) strncmp(x, mouse_id, length(mouse_id)), cellids_exp_valid)';
            
            if any(idx)
                % Get mean FRs of this mouse's neurons in this condition
                mouse_fr = expData(idx);
                
                % Compute stats
                mean_fr_mouse = mean(mouse_fr);
                se_fr_mouse = std(mouse_fr) / sqrt(length(mouse_fr));
                
                % Add to table
                newRow = array2table(...
                    {mouse_id, 'Experimental', condName, mean_fr_mouse, se_fr_mouse}, ...
                    'VariableNames', MouseTable.Properties.VariableNames);
                
                MouseTable = [MouseTable; newRow];
            end
        end
        
        % Process Control Mice
        for i = 1:length(mice_ctrl_list)
            mouse_id = mice_ctrl_list{i};
            
            % Find indices where cellid starts with mouse_id
            idx = cellfun(@(x) strncmp(x, mouse_id, length(mouse_id)), cellids_ctrl_valid)';
            
            if any(idx)
                % Get mean FRs of this mouse's neurons in this condition
                mouse_fr = ctrlData(idx);
                
                % Compute stats
                mean_fr_mouse = mean(mouse_fr);
                se_fr_mouse = std(mouse_fr) / sqrt(length(mouse_fr));
                
                % Add to table
                newRow = array2table(...
                    {mouse_id, 'Control', condName, mean_fr_mouse, se_fr_mouse}, ...
                    'VariableNames', MouseTable.Properties.VariableNames);
                
                MouseTable = [MouseTable; newRow];
            end
        end
        
    end
    
    % Initialize results table with appropriate column names
    ResultsTable = table([], [], [], [], [], [], 'VariableNames', ...
    {'Period', 'Mean_Exp', 'SE_Exp', 'Mean_Ctrl', 'SE_Ctrl', 'p value'}, ...
    'RowNames', {}); % Empty initially
    
    % Preallocate cell array for storing axes handles for later figure assembly
    axes_handles = cell(size(conditions, 1), 1);
    
    % Significance threshold
    alpha = 0.05;

    % Mean FR for conditions: bar plots and CDFs
    for cond = 1:size(conditions, 1)
        period_name = conditions{cond, 1}; % Name of condition
        d1 = conditions{cond, 2}; % Experimental data
        d2 = conditions{cond, 3}; % Control data
    
        % Plot bar graph and compute non-parametric statistical test (Mann-Whitney U); return means, standard errors, and p-value
        [mean_d1, se_d1, mean_d2, se_d2, Wp] = barmeanstat(d1, d2, WM_label, Ctrl_label, alpha, 'nonpaired', 0);
        
        current_condition = conditions{cond, 1};
        
        hold on;
        title(current_condition);
        ylabel('Mean FR (Hz)');
        hold off;
        axes_handles{cond*2 - 1} = gca; % Save bar plot axis handle
    
        % Plot CDF comparison between experimental and control groups
        cdf_fig(d1, d2, 100, WM_label, Ctrl_label, period_name);
        hold on;
        ylabel('Cumulative probability');
        xlabel('Mean FR (Hz)')
        hold off;
        axes_handles{cond*2} = gca; % Save CDF axis handle
    
        % Create a table row with the current condition's stats
        NewRow = table({period_name}, mean_d1, se_d1, mean_d2, se_d2, Wp, ...
            'VariableNames', ResultsTable.Properties.VariableNames);
    
        % Append to results table
        ResultsTable = [ResultsTable; NewRow];
    end
    
    % Display statistical summary in command window
    disp('Statistical Results:');
    disp(ResultsTable);
    
    % Save statistical results to Excel file
    filename = fullfile(resdir, 'StatResults.xlsx');
    sheetName = 'Summary';
    writetable(ResultsTable, filename, 'Sheet', sheetName, 'WriteMode', 'overwrite');
    
    % Get number of rows in ResultsTable + 2 for spacing
    startRow = height(ResultsTable) + 3; % +2 for spacing, +1 because Excel is 1-indexed
    
    % Add header for MouseBaselineTable
    mouseHeader = MouseTable.Properties.VariableNames;
    writecell(mouseHeader, filename, 'Sheet', sheetName, 'Range', ['A' num2str(startRow)]);
    
    % Write MouseTable starting just below ResultsTable
    startRow = startRow + 1;
    writecell(table2cell(MouseTable), filename, ...
             'Sheet', sheetName, 'Range', ['A' num2str(startRow)]);
    
    % Assemble composite figure (Fig 2E-J)
    
    % Define figure size in centimeters
    fig_width = 20;      
    fig_height = 10;   
    
    % Create figure and tiled layout
    figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height]);
    t = tiledlayout(2, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    existing_axes = axes_handles; % Axes handles from previous plots
    
    % Define normalized tick lengths for visual consistency
    tick_length_cm = 0.2;
    tick_length_y = tick_length_cm / fig_width;
    tick_length_x = tick_length_cm / fig_height;
    
    % Copy plots into tiled layout
    for iT = 1:12
        new_axes = nexttile(t, iT); % Create new tile
        children = get(existing_axes{iT}, 'Children'); % Get content of original axes
        copyobj(children, new_axes); % Copy content into new tile
    
        % Set visual styles
        set(new_axes, 'TickDir', 'out', 'Box', 'off', 'TickLength', [tick_length_x, tick_length_y]);
    
        % Preserve original axes labels
        ylabel(new_axes, get(get(existing_axes{iT}, 'YLabel'), 'String'));
        xlabel(new_axes, get(get(existing_axes{iT}, 'XLabel'), 'String'));
    
        if mod(iT,2) == 1
            new_axes.YLim = [0, 22]; % Standardize y-axis for bar plots
            new_axes.XTick = [1 2];
            new_axes.XTickLabel = {WM_label, Ctrl_label}; % Define tick labels
            new_axes.YTick = [0, 20]; % Show only 0 and 20 ticks
            new_axes.YLabel.VerticalAlignment = 'middle';
        else
            new_axes.YTick = [0, 1]; % For CDF plots, show only 0 and 1 ticks
        end
    
        % Add legend to CDF plots
        if mod(iT,2) == 0
            h = findobj(new_axes, '-regexp', 'DisplayName', '.+'); % Find legends in CDF original plots
            legend(new_axes, flipud(h), 'Location', 'southwest', 'Box', 'off'); % Add them to composite figure
        end
    end
    
    % Set font and line styles across all elements
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 9, 'LineWidth', 1);
    
    % Save final figure in both JPG and SVG formats
    fnm = [resdir '\' 'Stat panel - fig2.jpg'];
    saveas(gcf, fnm)
    
    fnmm = [resdir '\' 'Stat panel - fig2.svg'];
    saveas(gcf, fnmm)

end

function [mn_rew, mn_punish] = process_feedback_data(feedback_meanFR, normPSTH, type)
    % Filters feedback FR data based on valid rows in normPSTH, and returns
    % mean FRs for reward and punishment events

    % Identify rows that are fully NaN in normalized PSTH
    nan_rows = all(isnan(normPSTH), 2);  
    
    % Select correct subset from feedback FR data based on type
    if strcmp('Exp', type)
        mn_rew = feedback_meanFR{1,1}(~nan_rows, :); 
        mn_punish = feedback_meanFR{1,2}(~nan_rows, :); 
    else 
        mn_rew = feedback_meanFR{2,1}(~nan_rows, :); 
        mn_punish = feedback_meanFR{2,2}(~nan_rows, :); 
    end
end
