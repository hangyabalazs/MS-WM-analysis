function [groupCountsExp, groupCountsCtrl] = figure3_panel(datadir, resdir, cleaned_data, issave, varargin )
% FIGURE3_PANEL Generate figure panels comparing experimental and control neuron responses.
%   [GROUPCOUNTSEXP, GROUPCOUNTSCTRL] = FIGURE3_PANEL(DATADIR, RESDIR, CLEANED_DATA, ISSAVE, VARARGIN) 
%   generates Figure 3 comparing experimental and control neuron responses.
%   Chi-square test results are saved as an Excel file in the specified results directory.
%
%   Inputs:
%       datadir       - Directory containing input data files.
%       resdir        - Directory where results (figures, Excel files) will be saved.
%       cleaned_data  - Preprocessed data structure with PSTH and categorization matrices.
%       issave        - Logical flag to save figures (true/false).
%
%   Optional parameter-value pairs:
%       'delay'       - Delay period time window [start, end]. Default: [0.2, 1].
%       'window'      - Time window for analysis [start, end]. Default: [-2, 6].
%       'dt'          - Time resolution of bin raster in seconds. Default: 0.001.
%
%   Example:
%       figure3_panel('data/', 'results/',  cleaned_data, true);
%
%   Note:
%       Make sure the speaker image used in the figure is present in datadir.
%
% See also: GROUPING_PANEL_EXP, GROUPING_PANEL_CTRL.
%
%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025

    % Default arguments setup using inputParser
    prs = inputParser;
    addRequired(prs,'datadir',@(s)ischar(s));          % Directory where input data is stored
    addRequired(prs,'resdir',@(s)ischar(s));           % Directory to save resulting panels
    addRequired(prs,'cleaned_data',@(s)isstruct(s));   % Cleaned data structure (see process_ms_wm_data.m)
    addRequired(prs, 'issave');                         % Logical flag to save figures or not
    addParameter(prs, 'delay',[0.2 1]);                 % Delay period window [start end]
    addParameter(prs, 'window',[-2 6]);                  % Time window for PSTH data aligned to FixationBeginning
    addParameter(prs,'dt',0.001,@isnumeric);            % Time resolution of bin raster in seconds
    parse(prs,datadir,resdir,cleaned_data, issave,varargin{:});
    g = prs.Results;
    
    close all;                                          % Close all existing figures
    
    if ~isfolder(resdir)
        mkdir(resdir);                                  % Create results directory if it doesn't exist
    end
    
    % Extract normalized PSTH data for control and experimental groups
    normPSTH_ctrl = cleaned_data.spsth_data.normPSTH_ctrl_cl;
    normPSTH_exp = cleaned_data.spsth_data.normPSTH_exp_cl;
    
    % Define grouping matrices and labels 
    grouping_matrixCtrl = cleaned_data.delay_response.ResponseCategoriCtrl;
    grouping_matrixExp= cleaned_data.delay_response.ResponseCategoriExp;
    group_names = {'Inh','Act','Inh-Act','Act-Inh', 'NonResp'};
    labels = {'Inhibited','Activated','Inhibited then activated','Activated then inhibited','NonResponsive'};
    
    
    numGroups = 5;                                     % Number of groups/categories
    time = (g.window(1):g.dt:g.window(2));             % Time vector for PSTH
    
    try
        speaker_image = imread([datadir '\Speaker_picture.png']); % Load speaker image for figure overlay
    catch
        disp('Speaker image not found in data directory. Please make sure it is present before proceeding');
    end
    
    colors = {"#66CDAA","#EA9782","#0072BD", "#FF8383","#6161a7"}; % Colors for groups
    
    % Sort and compute statistics for experimental group PSTHs
    [averageExp, SEExp, groupCountsExp, psthsExp, groupboundariesExp, ~ ] = ...
        group_sort_psths(g, numGroups, group_names, grouping_matrixExp, normPSTH_exp, time);
    
    % Sort and compute statistics for control group PSTHs
    [averageCtrl, SECtrl, groupCountsCtrl, psthsCtrl, groupboundariesCtrl, time_cut] = ...
        group_sort_psths(g, numGroups, group_names, grouping_matrixCtrl, normPSTH_ctrl, time);
    
    filename = fullfile(datadir, 'StatResults.xlsx'); % Define Excel file path for saving stats
    
    emptyCol = repmat({''}, length(labels), 1);      % Empty column for formatting
    
    % Combine group labels and counts for experimental and control groups
    groupData = [cellstr(labels(:)), num2cell(groupCountsExp'), num2cell(groupCountsCtrl'), emptyCol];
    
    % Perform Chi-square tests 
    % Test: Exp neurons with significant change (groups 1-4) vs non-responsive (group 5)
    observed = [sum(groupCountsExp(1:4)), groupCountsExp(5)];
    expected = [sum(groupCountsExp)*0.05, sum(groupCountsExp)*0.95]; % Chance level expected counts
    [~, pval_chi, stats] = chi2gof(observed, 'Expected', expected, 'Emin', 5);
    fprintf('Chi-square goodness of fit test, Exp neurons significant change vs non-responsive\n');
    fprintf('Chi-square statistic: %.3f\n', stats.chi2stat);
    fprintf('Degrees of freedom: %.3f\n', stats.df);
    fprintf('Chi-square test p-value: %.4e\n', pval_chi);
    
    % Test: Exp vs Ctrl by primary response (3 categories)
    fprintf('Chi-square test, Exp vs Ctrl, grouped by primary response (3 categories)\n');
    [chi_stat_1, p_value_1, df_1] = chiSquareTest([groupCountsExp(1)+groupCountsExp(3), ...
        groupCountsExp(2)+groupCountsExp(4), groupCountsExp(5); ...
        groupCountsCtrl(1)+groupCountsCtrl(3), ...
        groupCountsCtrl(2)+groupCountsCtrl(4), groupCountsCtrl(5)], 0.001);
    
    % Test: Exp vs Ctrl full category distribution (5 categories)
    fprintf('Chi-square test, Exp vs Ctrl, (5 categories)\n');
    [chi_stat_2, p_value_2, df_2] = chiSquareTest([groupCountsExp; groupCountsCtrl], 0.001);
    
    % Prepare table for writing to Excel
    testRows = {
        'Stat grouping', 'Chi-square Statistic', 'Degrees of Freedom', 'p-value';
        'Exp neurons significant change vs non-responsive', stats.chi2stat, stats.df, pval_chi;
        'Exp vs Ctrl, grouped by primary response (3 categories)', chi_stat_1, df_1, p_value_1;
        'Exp vs Ctrl, full category distribution (5 categories)', chi_stat_2, df_2, p_value_2;
    };
    
    % Combine group counts and test results into one cell array
    combined = [
        {'Group Counts', '', '', ''};
        {'Group', 'Exp', 'Ctrl', ''};
        groupData;
        {'', '', '', ''};
        {'Chi-Square Test Results', '', '', ''};
        testRows
    ];
    
    % Write combined data to second sheet of Excel file
    writecell(combined, filename, 'Sheet', 2);
    
    
    % Create figure and plot panels
    figure;
    
    % Plot heatmap for experimental group
    subplot('Position',[0.09, 0.5838, 0.3, 0.3412]);
    plot_heatmap(time_cut, psthsExp, groupboundariesExp, 28, speaker_image, 0.141, 0.91);
    
    % Plot average PSTH for experimental group
    subplot('Position',[0.45, 0.5838, 0.2, 0.3412]);
    plot_avg(time_cut, averageExp, SEExp, colors, speaker_image, 0.501, 0.91);
    
    % Plot pie chart of group counts for experimental group
    subplot('Position',[0.72, 0.5838, 0.2074, 0.3412]);
    plot_piechart(groupCountsExp, 'Experimental', colors, labels);
    
    % Plot heatmap for control group
    subplot('Position',[0.09, 0.11, 0.3, 0.3412]);
    plot_heatmap(time_cut, psthsCtrl, groupboundariesCtrl, 28, speaker_image, 0.141, 0.4342);
    
    % Plot average PSTH for control group
    subplot('Position',[0.45, 0.11, 0.2, 0.3412]);
    plot_avg(time_cut, averageCtrl, SECtrl, colors, speaker_image, 0.501, 0.4362);
    
    % Plot pie chart of group counts for control group
    subplot('Position',[0.72, 0.03, 0.2074, 0.3412]);
    plot_piechart(groupCountsCtrl, 'Control', colors, labels);
    
    set(gcf, 'Renderer', 'painters'); % Use painters renderer for figure output
    
    % Define filename
    fnm = 'Fig3AB';
    
    % Save figures if requested
    if issave
        saveas(gcf, [resdir '\' fnm '_Stat.svg']);
        saveas(gcf, [resdir '\' fnm '_Stat.jpg']);
    end
    
    % Generate additional panels: Fig3C & S5C 
    grouping_panel_exp(resdir, datadir, grouping_matrixExp, normPSTH_exp);
    grouping_panel_ctrl(resdir, datadir, grouping_matrixCtrl, normPSTH_ctrl);

end

function [averages, SE, groupCounts, psths, groupboundaries, time_cut] = group_sort_psths(g,...
    numGroups, group_names, grouping_matrix, psth_data, time)
% Groups and sorts PSTHs based on neuron classifications,
% computes average PSTH and standard error (SE) for each group,
% and trims the PSTHs and time window for plotting.

    psths = cell(1, numGroups);
    averages = cell(1, numGroups);
    SE = cell(1, numGroups);
    groupCounts = zeros(1, numGroups);
    
    % Loop through each group to find member cells and calculate stats
    for group = 1:length(group_names)
        
        if isa(grouping_matrix,'double')
            % Numeric grouping: find indices directly
            groupIndices = find(grouping_matrix == group);
            groupCounts (group) = sum(grouping_matrix == group);
        else
            % Cell array of strings grouping: match group names
            groupLabel = group_names{group};
            groupIndices = find(strcmp(grouping_matrix, groupLabel));
            groupCounts = arrayfun(@(group) sum(strcmp(grouping_matrix,group)), group_names); 
        end
        
        % Extract PSTHs for group cells
        psths{group} = psth_data(groupIndices,:);
        % Compute average PSTH across cells in group
        averages{group} = mean(psths{group}, 1);
        % Compute standard error of the mean PSTH
        SE{group} = std(psths{group})/sqrt(size(psths{group}, 1));
    end

    sorted_psths = cell(1, numGroups);
    
    % Sort cells within each group by peak activity during delay period
    for group = 1:numGroups
        numCells = size(psths{group}, 1);
        mx = zeros(numCells, 1);
    
        for iCell = 1:numCells
            delay_start =g.delay(1);
            delay_end=g.delay(2);
            delay_idx = (time >= delay_start) & (time <= delay_end);
            % Find max firing rate within delay window for each cell
            mx(iCell) = max(psths{group}(iCell, delay_idx), [], 2);
        end
    
        % Sort indices by peak activity during the delay period
        [~, srtinx] = sort(mx, 'descend');
        sorted_psths{group} = psths{group}(srtinx, :);
    end
    
    % Concatenate sorted PSTHs from all groups into one matrix
    psths = nan(length(grouping_matrix), size(time, 2));
    psths(1:groupCounts(1), :) = sorted_psths{1};
    psths(groupCounts(1)+1:sum(groupCounts(1:2)), :) = sorted_psths{2};
    psths(sum(groupCounts(1:2))+1:sum(groupCounts(1:3)), :) = sorted_psths{3};
    psths(sum(groupCounts(1:3))+1:sum(groupCounts(1:4)), :) = sorted_psths{4};
    psths(sum(groupCounts(1:4))+1:end, :) = sorted_psths{5};
    
    % Calculate group boundary indices for plotting separation lines
    groupboundaries = cumsum(groupCounts);
    groupboundaries = groupboundaries(1:end-1);
    
    % Trim time window & psth matrix to focus on relevant period (-0.504 to 1.505 sec)
    [~, time_start_index] = min(abs(time - (-0.504)));
    [~, time_end_index] = min(abs(time - 1.505));
    time_cut = time(:,time_start_index:time_end_index);
    psths = psths(:, time_start_index:time_end_index);
    
    % Trim averages and SE matrices accordingly
    for iG = 1:numGroups
        SE{iG} = SE{iG}(time_start_index:time_end_index);
        averages{iG} = averages{iG}(time_start_index:time_end_index);
    end

end


function plot_heatmap(time, psths, groupboundaries, width, speaker_image, speaker_x, speaker_y)
% Plots heatmap of sorted PSTHs with group separation lines and overlays a speaker image.

    imagesc(time, 1:size(psths, 1), psths);
    ylabel('Neuron #');
    ylim([-width, size(psths, 1) + 0.5]);
    colorbar;
    clim([-9.51373008123223 9.51373008123223]); % Set color axis limits
    c = colorbar;
    c.Label.String = 'Normalized FR';
    c.Label.Rotation = -90;
    c.TickDirection = 'out';
    xlabel('Time from cue onset (s)');
    hold on;
    
    % Draw horizontal lines to separate groups
    for iC = 1:length(groupboundaries)
        line([time(1), time(end)], [groupboundaries(iC), groupboundaries(iC)], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    end
    
    % Draw fading rectangles and speaker image overlay
    plot_fading_rectangles([1, size(psths, 1)], speaker_image, -width, width, speaker_x, speaker_y);
    
    hold off;
end 

function plot_avg(time, average, SE, colors, speaker_image,  speaker_x,  speaker_y)
% Plots average PSTH with shaded standard error for each neuron group.

    hold on
    for iC = 1:length(average)
        errorshade(time, average{iC}, SE{iC}, 'LineColor', colors{iC}, 'ShadeColor', 'black', 'FaceAlpha', 0.4);
    end
    ylim([-6, 6]);
    xlim([time(1), time(end)]);
    xlabel('Time from cue onset (s)');
    ylabel('Average SPSTH');
    
    % Title depending on speaker image vertical position
    if speaker_y == 0.91
        title('DR-2AFC WM');
    else 
        title('DR-2AFC Control');
    end
    
    % Draw fading rectangles and speaker image overlay
    plot_fading_rectangles([-6, 5.5], speaker_image, 5.5, 0.5, speaker_x, speaker_y);
    
    hold off;
end

function plot_fading_rectangles(y_limits, speaker_image, y_pos, height, speaker_x, speaker_y)
% Draws a colored bar with fading transparency to highlight stimulus period
% and overlays a speaker image on the plot.

    fade_length = 50; % Length of fading effect
    transparency_gradient = linspace(0.4, 0.02, fade_length);
    x_pos=linspace(1,1.5,fade_length);
    start_index = 0.2;
    end_index = 1;
    
    % Base colors for solid and fading parts (RGBA for fading)
    solid_color=[0.811764705882353, 0.250980392156863, 0.513725490196078];
    fade_color = [0.811764705882353, 0.250980392156863, 0.513725490196078, 0.5]; 
    
    % Plot vertical lines for cue onset and markers
    line([0, 0], [y_limits(1), y_limits(2)], 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1);
    if speaker_y == 0.91
        line([0.2, 0.2], [y_limits(1), y_limits(2)], 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1);
    end
        
    % Plot solid colored rectangle bar representing stimulus period
    rectangle('Position', [start_index, y_pos, end_index - start_index, height], ...
        'FaceColor', solid_color, 'EdgeColor', solid_color);
    
    % Plot fading effect with transparency gradient
    for iC = 1:fade_length
        fade_color(4) = transparency_gradient(iC);  % Update transparency
        rectangle('Position', [x_pos(iC), y_pos, 0.01, height], 'FaceColor', fade_color, 'EdgeColor', fade_color);
    end
        
    % Format axes
    ax = gca;
    ax.TickDir = 'out';
    ax.Box = 'off';
        
    % Insert speaker image in the specified position
    ax2 = axes('Position', [speaker_x, speaker_y, 0.02, 0.02]);  
    imshow(speaker_image);
    axis('off');
    grid on;

end

function plot_piechart(groupCounts, dataset, colors, labels)
% Creates a pie chart to show proportions of neuron response categories.

    numGroups=size(groupCounts,2);
    explode = ones(1, numGroups); % Explode all slices equally
    
    h = pie(groupCounts, explode);
    
    % Color each pie slice according to group colors
     for iC = 1:size(groupCounts,2)
         h(2 * iC - 1).FaceColor = colors{iC};
     end
    
     % Add legend for Control dataset on specified position
     if strcmp(dataset,'Control')
        legend(labels, 'Position', [0.700691944490707,0.414854648652314,0.269642851821013,0.169047614435355], 'Box', 'off');
     end
    
    hold off;

end