function theta_piecharts_supp(resp_categ_ctrl, resp_categ_exp, theta_inx_ctrl,...
    theta_inx_exp, normPSTH_ctrl, normPSTH_exp, time, resdir)
%THETA_PIECHARTS_SUPP Plots S8A-B and theta-group-based pie charts (S9).
%   THETA_PIECHARTS_SUPP(RESP_CATEG_CTRL, RESP_CATEG_EXP, THETA_INX_CTRL,...
%   THETA_INX_EXP, NORMPSTH_EXP, NORMPSTH_CTRL, TIME, RESDIR)
%   plots distribution of firing rates and theta index and the proportion of
%   strongly, moderately, and non-theta-rhythmic neurons across response 
%   categories for both control and experimental datasets.
%
%   INPUTS:
%       resp_categ_ctrl - Cell array of response categories (control)
%       resp_categ_exp  - Cell array of response categories (experimental)
%       theta_inx_ctrl  - Theta index values for control neurons
%       theta_inx_exp   - Theta index values for WM neurons
%       normPSTH_ctrl   - Normalized PSTH matrix for control neurons
%       normPSTH_exp    - Normalized PSTH matrix for WM neurons
%       resdir          - Directory path to save output figures
%
%   See also: ACGMOD, CHISQUARETEST

%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025

    % Define group labels and plot appearance
    group_names = {'Inhibited','Activated','Non-responsive'}; % response categories' names
    colors_groups = {"#66CDAA","#EA9782","#6161a7"}; % colors for reponse categories
    colors_theta = {"#EBDB34","#3495EB","#A07BE0"}; % colors for theta/non-theta groups
    colors3s = [0.9215686274509803 0.8588235294117647 0.20392156862745098;
    0.20392156862745098 0.5843137254901961 0.9215686274509803;
    0.6274509803921569 0.4823529411764706 0.8784313725490196]; % colors for theta/non-theta groups in RGB
    explode3 = ones(1,3); % for piecharts
    theta_threshold_high = 0.5; % threshold for strong theta rhythmicity
    theta_threshold_mod = 0.2; % threshold for moderate theta rhythmicity

    % Process Experimental and Control datasets
    [~, theta_cells_exp, theta_proportions_exp, high_inx_exp, mod_inx_exp,...
        low_inx_exp] = process_theta_groups(resp_categ_exp, theta_inx_exp, ...
        theta_threshold_high, theta_threshold_mod);
    [~, theta_cells_ctrl, theta_proportions_ctrl, high_inx_ctrl, mod_inx_ctrl,...
        low_inx_ctrl] = process_theta_groups(resp_categ_ctrl, theta_inx_ctrl,...
        theta_threshold_high, theta_threshold_mod);

    % Define temporal parameters for delay
    delay_start = 0.2; 
    delay_end = 1; 
    
    % Compute indices for different temporal windows
    delay_inx = find(time >= delay_start & time <= delay_end);

    % Get mean FR during delay for all cells
    delay_mean_FR_exp = mean(normPSTH_exp(:, delay_inx), 2);
    delay_mean_FR_ctrl = mean(normPSTH_ctrl(:, delay_inx), 2);

    % Scatter plot Theta Index vs Mean FR (colored by group)
    figure('Units','centimeters', 'Position', [0,0, 20.45, 26.445]);
    tiledlayout(6, 4, 'TileSpacing', 'compact', 'Padding', 'compact');   
    nexttile([2 2]);
    hold on;
    scatter([theta_inx_exp(high_inx_exp); theta_inx_ctrl(high_inx_ctrl)], [delay_mean_FR_exp(high_inx_exp); delay_mean_FR_ctrl(high_inx_ctrl)], ...
        10, colors3s(1,:), 'filled');
    scatter([theta_inx_exp(mod_inx_exp); theta_inx_ctrl(mod_inx_ctrl)], [delay_mean_FR_exp(mod_inx_exp); delay_mean_FR_ctrl(mod_inx_ctrl)], ...
        10, colors3s(2,:), 'filled');
    scatter([theta_inx_exp(low_inx_exp); theta_inx_ctrl(low_inx_ctrl)], [delay_mean_FR_exp(low_inx_exp); delay_mean_FR_ctrl(low_inx_ctrl)], ...
        10, colors3s(3,:), 'filled');
    xlabel('Theta Index');
    ylabel('Mean FR - Delay (Hz)');
    legendEntries = {
        'Strongly theta-rhythmic', 'Moderately theta-rhythmic', 'Non-theta-rhythmic'};
    legend(legendEntries, 'Location', 'best', 'Box', 'off');
    xlim([-0.6 1]);
    ylim([-6 12]);
    set(gca, 'TickDir', 'out', 'Box', 'off');

    % Histogram of Theta Index 
    nexttile([2 2]);
    hold on;
    edges = linspace(0, 1, 21); % 20 bins
    % Experimental and ctrl, pooled
    h1 = histogram([theta_inx_exp(high_inx_exp) ; theta_inx_ctrl(high_inx_ctrl)], edges, 'FaceColor', colors_theta{1});
    h2 = histogram([theta_inx_exp(mod_inx_exp); theta_inx_ctrl(mod_inx_ctrl)] , edges, 'FaceColor', colors_theta{2});
    h3 = histogram([theta_inx_exp(low_inx_exp); theta_inx_ctrl(low_inx_ctrl)], edges, 'FaceColor', colors_theta{3});
    xlabel('Theta Index');
    ylabel('Number of Neurons');   
    legend([h1, h2, h3], legendEntries, ...
           'Location', 'best', 'Box', 'off');
    set(gca, 'TickDir', 'out', 'Box', 'off');

    % Run chi-square test on group-level theta counts
    chiSquareTest([theta_cells_exp; theta_cells_ctrl], 0.05); 

    % Plot the grouped pie charts
    plot_theta_piecharts(theta_cells_exp, theta_cells_ctrl,theta_proportions_exp,...
        theta_proportions_ctrl,group_names, colors_groups, colors_theta, explode3);

    % Save figure
    set(gcf, 'Renderer', 'painters');
    saveas(gcf, fullfile(resdir, 'Theta groups - S7.jpg'));
    saveas(gcf, fullfile(resdir, 'Theta groups - S7.svg'));
end

function [ratios, group_counts, theta_proportions, high_inx, mod_inx,...
    low_inx] = process_theta_groups(groups, theta_inx, high_thresh, mod_thresh)
% Computes theta proportions within each response category

    % Categorize each neuron based on theta index
    [theta_group, high_inx, mod_inx ,low_inx] = categorize_theta(theta_inx, high_thresh, mod_thresh);

    % Logical index of all rhythmic neurons (strong + moderate)
    thetaRhythmic = theta_inx >= mod_thresh;

    % Merged into 3: Inh (Inh + Inh-Act), Act (Act + Act-Inh), NonResp
    % Map original group indices to merged groups
    merge_map = [1, 2, 1, 2, 3]; % old index to new index

    % If cellstr, map based on original group names
    orig_names = {'Inh','Act','Inh-Act','Act-Inh', 'NonResp'};
    new_grouping = cell(size(groups));
    for i = 1:length(groups)
        idx = find(strcmp(orig_names, groups{i}));
        new_group = merge_map(idx);
        switch new_group
            case 1, new_grouping{i} = 'Inh';
            case 2, new_grouping{i} = 'Act';
            case 3, new_grouping{i} = 'NonResp';
        end
    end
    groups = new_grouping;

    % Updated group names and labels for 3 groups
    group_names = {'Inh','Act','NonResp'};

    % Count total number of cells per group
    totalPerGroup = arrayfun(@(g) sum(strcmp(groups, g)), group_names);

    % Count number of theta-rhythmic neurons per group
    group_counts = arrayfun(@(g) sum(strcmp(groups(thetaRhythmic), g)), group_names);

    % Compute theta ratio per group
    ratios = group_counts ./ totalPerGroup;

    % Compute proportions of theta strength subgroups within each response category
    nThetaGroups = 3;
    theta_proportions = zeros(3, nThetaGroups);
    for g = 1:3
        idx = strcmp(groups, group_names{g});
        for t = 1:nThetaGroups
            theta_proportions(g, t) = sum(theta_group(idx) == t);
        end
        theta_proportions(g,:) = theta_proportions(g,:) / sum(theta_proportions(g,:));
    end
end

function [thetaGroup, highInx, midInx, lowInx] = categorize_theta(theta_index, high_thresh, mod_thresh)
% Categorizes neurons into theta strength groups

  % Find all neuron indices
    allIndices = 1:length(theta_index);

    thetaGroup = zeros(size(theta_index));
    thetaGroup(theta_index > high_thresh) = 1;                               % Strongly theta-rhythmic
    thetaGroup(theta_index <= high_thresh & theta_index > mod_thresh) = 2;        % Moderately theta-rhythmic
    thetaGroup(theta_index <= mod_thresh) = 3;                              % Non-theta rhythmic

     % Categorize neurons based on Theta Index thresholds
    highInx = allIndices(theta_index > high_thresh);
    midInx = allIndices(theta_index < high_thresh & theta_index > mod_thresh);
    lowInx = allIndices(theta_index < mod_thresh);
end

function plot_theta_piecharts(exp_counts, ctrl_counts, prop_exp, prop_ctrl, group_names, colors5, colors3, explode2)
% Generates tiled pie chart layout for theta analysis.
   

     % Big pies for total theta rhythmic cells 
    nexttile(9, [2 1]);  % Experimental - spans 2 rows, 1 column
    plot_colored_pie(exp_counts , explode2, colors5);
    title('DR-2AFC WM');
    legend(group_names, "Box", "off", 'Position',[0.006250001064369,0.479412181270961,0.162499997871263,0.169047614435355]);

    % Small pies for each group (Experimental - top row)
    for g = 1:3
        nexttile(g + 9);  % Tiles 2 to 6
        plot_colored_pie(prop_exp(g,:), explode2, colors3);
        title({group_names{g}; ''});
    end
    for i=14:16
        nexttile(i);
        axis off;
    end

    nexttile;
    axis off;

    nexttile(17, [2 1]);  % Control (positioned in 3rd row, 1st col)
    plot_colored_pie(ctrl_counts, explode2, colors5);
    title('DR-2AFC Control');
    legend(group_names, "Box", "off", 'Position',[0.00547619047619,-0.003273809637342,0.162499997871263,0.169047614435355]);

    % Small pies for each group (Control - bottom row)
    for g = 1:3
        nexttile(g + 17);  % Tiles 12 to 16
        plot_colored_pie(prop_ctrl(g,:), explode2, colors3);
        title({group_names{g}; ''});
    end

    for i=22:24
        nexttile(i);
        axis off;
    end
end

function plot_colored_pie(data, explode, colors)
% Draw a pie chart with custom colors.

    h = pie(data, explode);
    for i = 1:length(data)
        h(2 * i - 1).FaceColor = colors{i};
    end
end

