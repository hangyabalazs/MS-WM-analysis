function stat_panel_fig4(cleaned_data, ROC, ROCtime, pvalues, resdir, roc_analysis_type)
%STAT_PANEL_FIG4 Generates statistical comparison panels for neuronal firing rates.
%   STAT_PANEL_FIG4(cleaned_data, ROC, ROCtime, pvalues, resdir) processes and
%   visualizes firing rate data comparing experimental (Exp) and control (Ctrl)
%   groups across conditions. It produces bar plots, cumulative density function
%   (CDF) plots, and ROC plots analyzing statistical differences during the delay
%   period.
%
%   Inputs:
%       cleaned_data : Struct with processed spiking data and statistics.
%       ROC          : ROC data for analysis.
%       ROCtime      : Time points corresponding to ROC.
%       pvalues      : P-values from ROC statistical comparisons.
%       resdir       : Directory path to save output figures and files.
%
%   Outputs:
%       Saves plots (.fig, .svg) and statistical results (.xlsx) to resdir.
% 
% See also : ROC_ANALYSIS, BARMEANSTAT, CDF_FIG, MEANBARGROUPS.
%
%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025
  
    close all 
    
    % Load and extract relevant variables from input struct
    ResponseCategoriExp=cleaned_data.delay_response.ResponseCategoriExp;
    ResponseCategoriCtrl=cleaned_data.delay_response.ResponseCategoriCtrl;
    normPSTH_exp=cleaned_data.spsth_data.normPSTH_exp;
    normPSTH_ctrl=cleaned_data.spsth_data.normPSTH_ctrl;
    spsthExp_raw=cleaned_data.spsth_data.spsthExp_raw;
    spsthCtrl_raw=cleaned_data.spsth_data.spsthCtrl_raw;
    statsExp=cleaned_data.stats.stats_Exp_cl;
    statsCtrl=cleaned_data.stats.stats_Ctrl_cl;
    time=cleaned_data.time;
    spsthExp = spsthExp_raw(sum(isnan(normPSTH_exp),2)==0,:);
    spsthCtrl = spsthCtrl_raw(sum(isnan(normPSTH_ctrl),2)==0,:);
    
    % Statistical significance threshold
    alpha = 0.05;
    
    % Define labels and groups
    label1='WM';
    label2='Ctrl';
    %group_names2={'Inh','Act','Inh-Act','Act-Inh', 'Non-resp'};
    groupAll={'Inh Exp','Inh Ctrl','Act Exp','Act Ctrl','Inh-Act Exp','Inh-Act Ctrl','Act-Inh Exp','Act-Inh Ctrl', 'NonResp Exp','NonResp Ctrl'};
    
    % Convert response categories to categorical arrays
    groupLabelsCtrl = categorical(ResponseCategoriCtrl);
    groupLabelsExp = categorical(ResponseCategoriExp);
    
    % Process stats data for WM neurons
    [nor_maxExp, nor_minExp, norFR1sthalfExp, norFR2ndhalfExp, nor_max_groupsExp, nor_min_groupsExp, ...
        nor_spsth_groupsExp, maxExp, minExp, FR1sthalfExp, FR2ndhalfExp,...
        min_groupsExp, max_groupsExp] = process_stats(time, spsthExp, statsExp, groupLabelsExp);
    
    % Process stats data for control neurons
    [nor_maxCtrl, nor_minCtrl, norFR1sthalfCtrl, norFR2ndhalfCtrl, nor_max_groupsCtrl, nor_min_groupsCtrl, ...
        nor_spsth_groupsCtrl, maxCtrl, minCtrl, FR1sthalfCtrl, FR2ndhalfCtrl,...
        min_groupsCtrl, max_groupsCtrl] = process_stats(time, spsthCtrl, statsCtrl, groupLabelsCtrl);
    
    % Organize 1st half, 2nd half spsth data, min/max groups into one cell array for group comparison
    [nor_spsth1sthalf, nor_spsth2ndhalf, mingroups, maxgroups] = organize_data(time,...
        nor_spsth_groupsExp, nor_spsth_groupsCtrl, min_groupsExp, min_groupsCtrl, max_groupsExp, max_groupsCtrl);
    
    % Organize min/max groups normalized spsth for baseline-subtracted analysis
    [~, ~, nor_min_groups, nor_max_groups] = organize_data(time, nor_spsth_groupsExp, nor_spsth_groupsCtrl,...
        nor_min_groupsExp, nor_min_groupsCtrl, nor_max_groupsExp, nor_max_groupsCtrl);
    
    % Plotting & statistical testing
    % Plot 1: Bar graph Min FR Exp vs Ctrl 
    [M1, SE1, M2, SE2, Wp1] = barmeanstat(nor_minExp, nor_minCtrl, label1, label2, alpha,'nonpaired', 0);    
    ylabel('Mean min FR (Hz)');
    title('All neurons');
    axes1=gca;
    
    % Plot 2: CDF Min FR Exp vs Ctrl
    cdf_fig(nor_minExp, nor_minCtrl, 100, label1, label2, 'Min FR');
    ylabel('Cumulative P');
    xlabel('Mean min FR (Hz)')
    axes2=gca;
    axes2.Title.String='';
    
    % Plot 3: Bar graph Max FR Exp vs Ctr
    [M3, SE3, M4, SE4, Wp2] = barmeanstat(nor_maxExp, nor_maxCtrl, label1, label2, alpha,'nonpaired', 0);    
    ylabel('Mean max FR (Hz)');
    title('All neurons');
    axes3=gca;
    
    % Plot 4: CDF Max FR Exp vs Ctrl
    [axes4,~] = cdf_fig(nor_maxExp, nor_maxCtrl, 100, label1, label2, 'Max FR');
    ylabel('Cumulative P');
    xlabel('Mean max FR (Hz)')
    axes4.Title.String='';
    
    % Plot 5: Bar graph Min FR Inh Exp vs Inh Ctrl
    [MinInhE, SEMIE, MinInhC, SEMIC, WpMI] = barmeanstat(nor_min_groupsExp{1}, nor_min_groupsCtrl{1}, label1, label2, alpha,'nonpaired', 0); 
    ylabel('Mean min FR (Hz)');
    title('Inhibited neurons');
    axes5=gca;
    
    % Plot 6: Bar graph Max FR Act Exp vs Act Ctrl
    [MaxActE, SEMAE, MaxActC, SEMAC, WpMA] = barmeanstat(nor_max_groupsExp{2}, nor_max_groupsCtrl{2}, label1, label2, alpha,'nonpaired', 0);
    ylabel('Mean max FR (Hz)');
    title('Activated neurons');
    axes6=gca;
    
    % Plot 7: Bar graph Mean FR 1st half delay Exp vs Ctrl (all neurons)
    [M5, SE5, M6, SE6, Wp3] = barmeanstat(norFR1sthalfExp, norFR1sthalfCtrl, label1, label2, alpha,'nonpaired', 0);    
    title('1st Half - Delay');
    ylabel('Mean FR (Hz)');
    axes7=gca;
    
    % Plot 8: Bar graph Mean FR 2nd half delay Exp vs Ctrl (all neurons)
    [M7, SE7, M8, SE8, Wp4] = barmeanstat(norFR2ndhalfExp, norFR2ndhalfCtrl, label1, label2, alpha,'nonpaired', 0);  
    title('2nd Half - Delay');
    ylabel('Mean FR (Hz)');
    axes8=gca;
    
    % Plot 9: Bar graph Mean FR 1st half delay per category
    [axes9, ~, Wp5, means1, SEs1] = meanbargroups(nor_spsth1sthalf, groupAll, alpha, '1st Half - Delay', 0);
    ylabel('Mean FR (Hz)');
    
    % Plot 10: Bar graph Mean FR 2nd half delay per category
    [axes10, ~, Wp6, means2, SEs2] = meanbargroups(nor_spsth2ndhalf, groupAll, alpha, '2nd Half - Delay', 0);
    ylabel('Mean FR (Hz)');
    
    % Create table for statistical results summary (means, SE, p-values)
    resultsTable1 = table(...
        M1, SE1, M2, SE2, Wp1, ...
        M3, SE3, M4, SE4, Wp2, ...
        MinInhE, SEMIE, MinInhC, SEMIC, WpMI, ...
        MaxActE, SEMAE, MaxActC, SEMAC, WpMA,...
        M5, SE5, M6, SE6, Wp3, ...
        M7, SE7, M8, SE8, Wp4, ...
        'VariableNames', { ...
        'MinFR_Exp', 'SE_MinFR_Exp', 'MinFR_Ctrl', 'SE_MinFR_Ctrl', 'Wp_MinFR', ...
        'MaxFR_Exp', 'SE_MaxFR_Exp', 'MaxFR_Ctrl', 'SE_MaxFR_Ctrl', 'Wp_MaxFR', ...
        'MinInhExp', 'SE_MinInhExp', 'MinInhCtrl', 'SE_MinInhCtrl', 'WpMI', ...
        'MaxActExp', 'SE_MaxActExp', 'MaxActCtrl', 'SE_MaActCtrl', 'WpMA', ...
        '1stHalfFR_Exp', 'SE_1stHalfFR_Exp', '1stHalfFR_Ctrl', 'SE_1stHalfFR_Ctrl', 'Wp_1stHalfFR', ...
        '2ndHalfFR_Exp', 'SE_2ndHalfFR_Exp', '2ndHalfFR_Ctrl', 'SE_2ndHalfFR_Ctrl', 'Wp_2ndHalfFR'...
        });
    
    resultsTable2 = table(...
        means1, SEs1, ...
        means2, SEs2, ...
        'VariableNames', { ...
        'Means_1stHalfGroups', 'SEs_1stHalfGroups', ...
        'Means_2ndHalfGroups', 'SEs_2ndHalfGroups' ...
        });
    
    % Save group p-value matrices separately
    % Make sure the Excel file is NOT open in another program before running writematrix, otherwise MATLAB will throw an error.
    writematrix(Wp5, [resdir '\FiringRateResults.xlsx'], 'Sheet', 'Wp_1stHalfGroups');
    writematrix(Wp6, [resdir '\FiringRateResults.xlsx'], 'Sheet', 'Wp_2ndHalfGroups');
    
    % Save tables to Excel
    writetable(resultsTable1, [resdir '\FiringRateResults.xlsx'], 'Sheet', 'GroupedResults1');
    writetable(resultsTable2, [resdir '\FiringRateResults.xlsx'], 'Sheet', 'GroupedResults2');
    
    % Build stat summary panel figure
    
    fig_width = 9;  
    fig_height = 20; 
    
    % Create figure with specified size in centimeters
    fig = figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height]);
    
    % Define a tiled layout with 6 rows 
    tiledlayout(6, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % List of axes handles from previously generated plots
    existing_axes = {axes1, axes2, axes3, axes4, axes5, axes6, axes7, axes8, axes9, axes10};
    
    % Loop through first 10 axes and copy their contents into the tiled layout
    for i = 1:10
        nexttile(i);                        % Move to the i-th tile in layout
        new_axes = gca;                    % Get current axes in tile
        children = get(existing_axes{i}, 'Children');  % Get all plot elements from old axes
        copyobj(children, new_axes);      % Copy plot elements to new axes
        title(new_axes, get(get(existing_axes{i}, 'Title'), 'String'));  % Copy title
        set(new_axes, 'TickDir', 'out', 'Box', 'off', 'FontSize', 7);    % Standardize appearance
        ylabel(new_axes, get(get(existing_axes{i}, 'YLabel'), 'String')); % Copy Y-axis label
    
        % Set X-axis ticks and labels for grouped bar plots
        if i == 1 || i == 3 || i == 5 || i == 6 || i == 7 || i == 8
            new_axes.XTick = [1 2];
            new_axes.XTickLabel = {label1, label2};
        else 
            xlabel(new_axes, get(get(existing_axes{i}, 'XLabel'), 'String')); % Copy X-axis label
        end
    
        % Add legend for CDF plot
        if i == 2 || i == 4 
            h = findobj(new_axes, '-regexp', 'DisplayName', '.+'); % Find legends in CDF original plots
            legend(new_axes, flipud(h), 'Location', 'southwest', 'Box', 'off'); % Add them to composite figure
        end
    
        % Manually adjust Y-limits for specific plots for visual clarity
        if i == 3
            ylim([0, 8]);
        elseif i == 1
            ylim([-10, 0]);
        elseif i == 7 || i == 8
            ylim([-3.5, 0]);
        elseif i == 9 || i == 10
            ylim([-7, 7]);
            xticks(new_axes, get(existing_axes{i}, 'XTick'));  % Copy X-ticks
            xticklabels(new_axes, get(existing_axes{i}, 'XTickLabel')); % Copy labels
        end
    end
    
    % Add ROC plots in the last row
    nexttile(11);
    % plot_roc(gca, ROCtime, pvalues, ROC, 1); % ROC for Inh groups
    plot_roc_log(gca, ROCtime, pvalues, ROC, 1); % ROC for Inh groups
    
    nexttile(12);
    % plot_roc(gca, ROCtime, pvalues, ROC, 2); % ROC for Act groups
    plot_roc_log(gca, ROCtime, pvalues, ROC, 2); % ROC for Act groups
    
    % Standardize font size across the entire figure
    set(findall(fig, '-property', 'FontSize'), 'FontSize', 8);
    
    % Save final figure to both JPG and SVG formats
    fnm = [resdir '\'  'Stat panel - fig4' roc_analysis_type '.jpg'];   % save
    saveas(gcf,fnm)
    
    fnmm = [resdir '\'  'Stat panel - fig4' roc_analysis_type '.svg'];   % save
    saveas(gcf,fnmm)
    
    % Supplementary Figure S6
    stat_categories(cleaned_data, resdir);
    
    % Supplementary Figure S7
    close all  % Close all open figures to start clean
    
    % Plot 1: Minimum firing rate (FR) during delay (All neurons)
    barmeanstat(minExp, minCtrl, label1, label2, alpha, 'nonpaired', 0);
    title('All neurons');
    ylabel('Mean min FR (Hz)')
    ax1 = gca;
    
    % Plot 2: CDF of minimum FR
    [ax2, ~] = cdf_fig(minExp, minCtrl, 100, label1, label2, 'Min FR');
    ylabel('Cumulative P');
    xlabel('Mean min FR (Hz)')
    ax2.Title.String = '';
    
    % Plot 3: Maximum FR during delay (All neurons)
    barmeanstat(maxExp, maxCtrl, label1, label2, alpha, 'nonpaired', 0);
    title('All neurons');
    ylabel('Mean max FR (Hz)')
    ax3 = gca;
    
    % Plot 4: CDF of maximum FR
    [ax4, ~] = cdf_fig(maxExp, maxCtrl, 100, label1, label2, 'Max FR');
    ylabel('Cumulative P');
    xlabel('Mean max FR (Hz)')
    ax4.Title.String = '';
    
    % Plot 5: Grouped minimum FR by category
    [ax5, ~] = meanbargroups(mingroups, groupAll, alpha, 'Min', 0);
    ylabel('Mean min FR (Hz)');
    ax5.Title.String = '';
    
    % Plot 6: Grouped maximum FR by category
    [ax6, ~] = meanbargroups(maxgroups, groupAll, alpha, 'Max', 0);
    ylabel('Mean max FR (Hz)');
    ax6.Title.String = '';
    
    % Plot 7: Grouped min FR (baseline subtracted)
    [ax7, ~] = meanbargroups(nor_min_groups, groupAll, alpha, 'Min', 0);
    ylabel('Mean min FR (Hz)');
    ax7.Title.String = '';
    
    % Plot 8: Grouped max FR (baseline subtracted)
    [ax8, ~] = meanbargroups(nor_max_groups, groupAll, alpha, 'Max', 0);
    ylabel('Mean max FR (Hz)');
    ax8.Title.String = '';
    
    % Plot 9: First half of delay period FR comparison
    barmeanstat(FR1sthalfExp, FR1sthalfCtrl, 'WM', 'Ctrl',alpha, 'nonpaired', 0);
    title('1st Half - Delay');
    ylabel('Mean FR (Hz)');
    ax9 = gca;
    
    % Plot 10: CDF for first half delay period
    cdf_fig(FR1sthalfExp, FR1sthalfCtrl, 100, label1, label2, 'First');
    ylabel('Cumulative P');
    xlabel('Mean FR (Hz)');
    ax10 = gca;
    ax10.Title.String = '';
    
    % Plot 11: Second half of delay period FR comparison
    barmeanstat(FR2ndhalfExp, FR2ndhalfCtrl, 'WM', 'Ctrl', alpha, 'nonpaired', 0);
    title('2nd Half - Delay');
    ylabel('Mean FR (Hz)');
    ax11 = gca;
    
    % Plot 12: CDF for second half delay period
    cdf_fig(FR2ndhalfExp, FR2ndhalfCtrl, 100, label1, label2, 'Second');
    ylabel('Cumulative P');
    xlabel('Mean FR (Hz)');
    ax12 = gca;
    ax12.Title.String = '';
    
    % Initialize figure for tiled layout to assemble all subplots
    fig_width = 9; 
    fig_height = 24; 
    fig = figure('Units', 'centimeters', 'Position', [0, 0, fig_width, fig_height]);
    tiledlayout(7, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % List of axes to copy into tiled layout
    existing_axes = {ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12};
    
    % Copy each individual plot into the tiled layout
    for i = 1:12
        nexttile(i);
        new_axes = gca;
        children = get(existing_axes{i}, 'Children');
        copyobj(children, new_axes);
        
        % Copy title and Y-axis label
        title(new_axes, get(get(existing_axes{i}, 'Title'), 'String'));
        ylabel(new_axes, get(get(existing_axes{i}, 'YLabel'), 'String'));
        xlabel(new_axes, get(get(existing_axes{i}, 'XLabel'), 'String'));
        
        % Set visual properties
        set(new_axes, 'TickDir', 'out', 'Box', 'off', 'FontSize', 7);
    
        % Apply custom X-ticks for bar plots
        if ~(i == 2 || i == 4 || i == 10 || i == 12)  % Skip CDF plots
            xticks(new_axes, get(existing_axes{i}, 'XTick'));
            xticklabels(new_axes, get(existing_axes{i}, 'XTickLabel'));
        end
        
        % Set specific X-axis labels
        if i == 1 || i == 3
            new_axes.XTick = [1 2];
            new_axes.XTickLabel = {label1, label2};
        end
    
        % Add legend to CDF plots
        if i == 2 || i == 4 || i == 10 || i == 12
            h = findobj(new_axes, '-regexp', 'DisplayName', '.+'); % Find legends in CDF original plots
            legend(new_axes, flipud(h), 'Location', 'southwest', 'Box', 'off'); % Add them to composite figure
        end
    
        % Manually define Y-limits for consistency across plots
        switch i
            case 1
                ylim([0, 15]);
            case 3
                ylim([0, 30]);
            case 5
                ylim([0, 16]);
            case 6
                ylim([0, 40]);
            case 7
                ylim([-12, 0]);
            case 8
                ylim([0, 20]);
            case {9, 11}
                ylim([0, 20]);
        end
    end
    
    % Add final row: ROC plots for additional group comparisons 
    nexttile(13);
    plot_roc_log(gca, ROCtime, pvalues, ROC, 3);  % ROC for InhAct groups
    
    nexttile(14);
    plot_roc_log(gca, ROCtime, pvalues, ROC, 4);  % ROC for ActInh groups
    
    % Standardize font size across all elements
    set(findall(fig, '-property', 'FontSize'), 'FontSize', 8);
    
    % Save figure in both JPG and SVG formats 
    saveas(gcf, [resdir '\'  'Stat panel - S7' roc_analysis_type '.jpg']);
    saveas(gcf, [resdir '\'  'Stat panel - S7' roc_analysis_type '.svg']);

end


function [nor_max, nor_min, norFR1sthalf, norFR2ndhalf, nor_max_groups, nor_min_groups, ...
          nor_spsth_groups, max, min, FR1sthalf, FR2ndhalf, min_groups, max_groups] ...
          = process_stats(time, spsth, stats, grouping_matrix)
% Computes baseline-subtracted and raw firing rate statistics (min, max, mean)
% across time windows and groups for neural data.

% Define delay period indices: full delay, 1st half, and 2nd half
delay_start = 0.2;
firsthalf_end = 0.6;
delay_end = 1;

firsthalf_inx = find(time >= delay_start & time <= firsthalf_end);
secondhalf_inx = find(time >= firsthalf_end & time <= delay_end);

% Define group names and number of groups
group_names = {'Inh', 'Act', 'Inh-Act', 'Act-Inh', 'NonResp'};
num_groups = numel(group_names);

% Compute baseline-subtracted PSTH, max and min FR, and raw max/min
[nor_spsth, nor_max, nor_min, max, min] = get_stats(spsth, stats);

% Group baseline-subtracted max/min FR by response category
nor_max_groups = categorize_data(num_groups, nor_max, grouping_matrix, group_names);
nor_min_groups = categorize_data(num_groups, nor_min, grouping_matrix, group_names);

% Mean FR during 1st and 2nd half of delay (raw)
FR1sthalf = mean(spsth(:, firsthalf_inx), 2);
FR2ndhalf = mean(spsth(:, secondhalf_inx), 2);

% Mean FR during 1st and 2nd half of delay (baseline-subtracted)
norFR1sthalf = mean(nor_spsth(:, firsthalf_inx), 2);
norFR2ndhalf = mean(nor_spsth(:, secondhalf_inx), 2);

% Group full PSTH and raw min/max values
nor_spsth_groups = categorize_data(num_groups, nor_spsth, grouping_matrix, group_names);
min_groups = categorize_data(num_groups, min, grouping_matrix, group_names);
max_groups = categorize_data(num_groups, max, grouping_matrix, group_names);

end


function [spsth1sthalf, spsth2ndhalf, min_data, max_data] = organize_data(time, spsth_groups1, ...
    spsth_groups2, min_groups1, min_groups2, max_groups1, max_groups2)
% Merges and processes grouped neural data from two datasets to compute min, 
% max, and mean firing rates during early and late delay periods.

% Define time window indices for the full delay period, first half, and second half
delay_start = 0.2;
firsthalf_end = 0.6;
delay_end = 1;

firsthalf_inx = find(time >= delay_start & time <= firsthalf_end);
secondhalf_inx = find(time >= firsthalf_end & time <= delay_end);

% Determine the number of groups
num_groups = size(spsth_groups1, 2);

% Preallocate cell arrays to combine categorized data from both datasets
spsth_data = cell(1, num_groups * 2);
min_data = cell(1, num_groups * 2);
max_data = cell(1, num_groups * 2);

% Interleave group data: dataset1 followed by dataset2 for each group
for group = 1:num_groups
    rowIndex1 = (group - 1) * 2 + 1; % Index for dataset1
    rowIndex2 = (group - 1) * 2 + 2; % Index for dataset2

    spsth_data{rowIndex1} = spsth_groups1{group};
    spsth_data{rowIndex2} = spsth_groups2{group};

    min_data{rowIndex1} = min_groups1{group};
    min_data{rowIndex2} = min_groups2{group};

    max_data{rowIndex1} = max_groups1{group};
    max_data{rowIndex2} = max_groups2{group};
end

% Compute mean FR for the first and second halves of the delay period
spsth1sthalf = cell(1, num_groups * 2);
spsth2ndhalf = cell(1, num_groups * 2);

for i = 1:length(spsth_data)
    % Average firing rate over the first half of delay
    spsth1sthalf{i} = mean(spsth_data{i}(:, firsthalf_inx), 2);
    
    % Average firing rate over the second half of delay
    spsth2ndhalf{i} = mean(spsth_data{i}(:, secondhalf_inx), 2);
end


end

function [nor_spsth, nor_max, nor_min, maxvalues, minvalues] = get_stats(spsth, stats)
% Extracts and normalizes neural response statistics from PSTH data.

    % Extract stats from structure fields
    baseline = cell2mat({stats.baseline}.');     % Mean baseline firing rate
    maxvalues = cell2mat({stats.maxvalue}.');    % Maximum firing rate during delay
    minvalues = cell2mat({stats.minvalue}.');    % Minimum firing rate during delay

    % Normalize max and min firing rates by subtracting the baseline
    nor_max = maxvalues - baseline;
    nor_min = minvalues - baseline;

    % Normalize the full PSTH matrix by subtracting each neuron's baseline
    nor_spsth = spsth - baseline;

end

function group_data = categorize_data(num_groups, allneurons_data, grouping_matrix, group_names)
% Groups neuron data based on specified group labels.

    group_data = cell(1, num_groups);
    for g = 1:num_groups
        group_data{g} = allneurons_data(grouping_matrix == group_names{g}, :);
    end

end

function plot_roc_log(ax, ROCtime, pvalues, ROC, group)
% Plots ROC and p-values using a log y-scale for p-values.

    yyaxis(ax, 'left');
    p = pvalues(group, :);
    p(p == 0) = 1e-10; % Avoid log(0)
    bar(ax, ROCtime, p, 'FaceAlpha', 0.5);
    set(ax, 'YScale', 'log');
    ylabel(ax, 'p value (log scale)');
    if group==3 || group ==4
        ylim([0.0001 10]);
        yticks([0.0001 1 10]);
    end

    hold(ax, 'on');
    yline(ax, 0.01, '--', 'LineWidth', 1);

    % Plot ROC values on the right y-axis
    yyaxis(ax, 'right');
    plot(ax, ROCtime, ROC(group, :));
    ylabel(ax, 'ROC', 'Rotation', -90, 'VerticalAlignment', 'cap');
    yticks(ylim);

    xlabel(ax, 'Time (s)');
    xticks([0.2 0.6 1]);
    set(ax, 'TickDir', 'out', 'Box', 'off');

    titles = {'Inhibited neurons', 'Activated neurons', ...
              'Inh-Act neurons', 'Act-Inh neurons'};
    title(ax, titles{group});
end
