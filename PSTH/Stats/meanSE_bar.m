function [H, Wp, means, SEs] = meanSE_bar(data, labels, alpha, stat_type, bar_title, edge_color)
%MEANSE_BAR  Bar plot of group means with standard errors and pairwise nonparametric tests
%   [H, Wp, means, SEs] = MEANSE_BAR(DATA, LABELS) creates a bar plot displaying the mean and
%   standard error for each group in DATA, labeled using LABELS. Pairwise nonparametric
%   statistical tests are performed and significant comparisons are annotated on the plot.
%
%   [___] = MEANSE_BAR(DATA, LABELS, ALPHA) specifies the significance level for hypothesis
%   testing (default: 0.01).
%
%   [___] = MEANSE_BAR(___, STAT_TYPE) specifies the type of test:
%       'nonpaired' — Mann-Whitney U test (default, for independent samples)
%       'paired'    — Wilcoxon signed-rank test (for matched or repeated measures)
%
%   [___] = MEANSE_BAR(___, BAR_TITLE) sets the title of the figure.
%
%   [___] = MEANSE_BAR(___, EDGE_COLOR) specifies the edge color for each bar as a cell array
%   of color specifications (e.g., {'r', 'b'}); defaults to black if omitted.
%
%   Inputs:
%     DATA          — Group data, specified as either:
%                     An N-by-G numeric matrix (N observations, G groups), or
%                     A 1-by-G cell array where each cell contains a numeric vector of
%                       observations for one group.
%     LABELS        — Group labels, specified as a 1-by-G cell array of character vectors or
%                     strings.
%     ALPHA         — Significance level (scalar in (0,1); default: 0.01).
%     STAT_TYPE     — Test type ('nonpaired' or 'paired'; default: 'nonpaired').
%     BAR_TITLE     — Plot title (character vector, string, or empty; default: '').
%     EDGE_COLOR    — Bar edge colors, specified as a 1-by-G cell array of color specs
%                     (e.g., {'k', [0.2 0.4 0.8]}); default: {'k', 'k', ..., 'k'}.
%
%   Outputs:
%     H             — Figure handle.
%     Wp            — Matrix of pairwise p-values.
%     means         — Vector of group means.
%     SEs           — Vector of standard errors.
%
%   Example:
%       data = {randn(20,1), randn(20,1) + 0.8};
%       labels = {'Control', 'Experimental'};
%       [H, Wp, means, SEs] = meanSE_bar(data, labels, 0.05, 'nonpaired', ...
%           'FR comparison', {'k', 'k'});
%
%   See also SIGNRANK, RANKSUM.

%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025

% Input argument check
narginchk(2, 6);
if nargin < 4
    stat_type = 'nonpaired';
end
if nargin < 3 || isempty(alpha)
    alpha = 0.01;   % default significance level
    edge_color = {'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k'};
    bar_title = '';
end

% Calculate means and SEs
numGroups = size(data, 2);
means=nan(numGroups,1);
SEs=nan(numGroups,1);
for i=1:numGroups
    means(i) = nanmean(data{i});
    SEs(i) = nanstd(data{i}) ./ sqrt(size(data{i}, 1));
end


% Bar plot
H = figure;
hold on;
% Define spacing between groups
spacing = 1.5; % Space between pairs of bars

% Plot bars with the specified arrangement (pairs of bars without space, then a gap)
barPositions = zeros(numGroups, 1);   % Store x-positions for bars
for i = 1:numGroups
    % Calculate position for each bar, alternating between no space for pairs and spacing between them
    groupIndex = ceil(i / 2);  % Group index (1, 2, 3, ...)
    position = (groupIndex - 1) * (2 + spacing) + (mod(i, 2) == 0); % Adjust positions for pairs
    bar(position, means(i), 'FaceColor', 'none', 'EdgeColor', edge_color{i}, 'LineWidth', 1);
    errorbar(position, means(i), SEs(i), 'k', 'LineStyle', 'none', 'CapSize', 0);
    barPositions(i) = position; % Store the x-position for the bar
end
title(bar_title); % Plot title

% Perform pairwise statistical tests
Wp = zeros(numGroups);
Wh = false(numGroups);
for i = 1:numGroups
    for j = i+1:numGroups
        switch stat_type
            case 'nonpaired'
                [Wp(i, j), Wh(i, j)] = ranksum(data{i}, data{j}, 'alpha', alpha);
            case 'paired'
                [Wp(i, j), Wh(i, j)] = signrank(data{i}, data{j}, 'alpha', alpha);
            otherwise
                error('boxstat:inputArg', 'Unsupported input argument.')
        end
    end
end

% Adjust y-axis limits
m1 = min(means - SEs);
m2 = max(means + SEs);
y_limits = [m1 - 2, m2 + 2];
ylim(y_limits);

ax = gca;
ax.TickDir = 'out';
ax.Box = 'off';
ax.XTick = 1:numGroups;
ax.XTickLabel = labels;

% Identify significant comparisons
significantPairs = [1,2; 3,4; 5,6; 7,8; 9,10];
numSignificant = 0;
filteredPairs = [];
for k = 1:size(significantPairs, 1)
    i = significantPairs(k, 1);
    j = significantPairs(k, 2);
    if i <= numGroups && j <= numGroups && Wh(i, j)
        filteredPairs = [filteredPairs; i, j, Wp(i, j)];
        numSignificant = numSignificant + 1;
    end
end

% Calculate new y-limits based on the number of significant comparisons
if numSignificant > 0
    extra_space = (numSignificant + 1) * 3; % Adjust this multiplier as needed for spacing
else
    extra_space = 1;
end
y_limits = [m1 - 2, m2 + extra_space];
ylim(ax, y_limits);

% Get current axis limits
y_limits = ylim(ax);
y_range = diff(y_limits);

% Define vertical line length as a fraction of y-range (e.g., 2%)
vertical_line_length = 0.01 * y_range;

% Horizontal offset for connecting lines (bar width dependent)
horizontal_offset = 0.1;
star_gap = 0.0001 * y_range;              % Distance between star and horizontal line

% Loop over significant pairs
for pair_idx = 1:numSignificant
    i = filteredPairs(pair_idx, 1);
    j = filteredPairs(pair_idx, 2);
    p_val = filteredPairs(pair_idx, 3);

     % Determine significance stars
    if p_val < 0.0001
        star = '***';
    elseif p_val < 0.001
        star = '**';
    elseif p_val < alpha
        star = '*';
    else
        continue;
    end

     % Position x between two bars
    tpos1 = mean([barPositions(i), barPositions(j)]);

    if means(i) >= 0 && means(j) >= 0
        % Bars are positive  draw lines above
        m1 = max(means(i)+SEs(i), means(j)+SEs(j));
        y_base = m1 + 0.005 * y_range; % margin above top bar

        % Top of horizontal line
        y_top = y_base;

        % Star position
        y_star = y_top + star_gap;

        % Vertical lines 
        y_bottom = y_top - vertical_line_length;

        % Draw horizontal line
        line([barPositions(i)+horizontal_offset, ...
                   barPositions(j)-horizontal_offset], ...
                  [y_top, y_top], 'Color', 'black', 'Parent', ax);

        % Left vertical line: up from line
        line([barPositions(i)+horizontal_offset, ...
                    barPositions(i)+horizontal_offset], ...
                   [y_top, y_bottom], 'Color', 'black', 'Parent', ax);

        % Right vertical line: up from line
        line([barPositions(j)-horizontal_offset, ...
                    barPositions(j)-horizontal_offset], ...
                   [y_top, y_bottom], 'Color', 'black', 'Parent', ax);

        % Add star above the horizontal line
        text(tpos1, y_star, star, ...
                   'HorizontalAlignment','center', ...
                   'VerticalAlignment','bottom', ...
                   'FontSize',12, 'Color','black', 'Parent', ax);

    else
        % Bars may be negative  draw line **below**
        m1 = min(means(i)-SEs(i), means(j)-SEs(j));
        y_base = m1 - 0.005 * y_range; % margin below lowest point

        % Bottom of horizontal line
        y_top = y_base;

        % Star
        y_star = y_top + star_gap;

        % Vertical 
        y_bottom = y_top + vertical_line_length;

        % Draw horizontal line
        line([barPositions(i)+horizontal_offset, ...
                   barPositions(j)-horizontal_offset], ...
                  [y_top, y_top], 'Color', 'black', 'Parent', ax);

        % Left vertical line: down from line
        line([barPositions(i)+horizontal_offset, ...
                    barPositions(i)+horizontal_offset], ...
                   [y_top, y_bottom], 'Color', 'black', 'Parent', ax);

        % Right vertical line: down from line
        line([barPositions(j)-horizontal_offset, ...
                    barPositions(j)-horizontal_offset], ...
                   [y_top, y_bottom], 'Color', 'black', 'Parent', ax);

        % Add star below the horizontal line
        text(tpos1, y_star, star, ...
                   'HorizontalAlignment','center', ...
                   'VerticalAlignment','top', ...
                   'FontSize',12, 'Color','black', 'Parent', ax);
    end
end
hold off;
end
