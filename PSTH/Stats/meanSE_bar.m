function [H, Wp, means, SEs] = meanSE_bar(data, labels, alpha, str, t, c)
% MEANSE_BARSTAT   Mean bar plot with statistics for multiple groups.
%
%   [H, WP] = MEANSE_BARSTAT3(DATA, LABELS, ALPHA, STR, T, C) plots mean bar plot (H, figure handle)
%   of data sets for multiple groups using labels. It performs pairwise
%   non-parametric tests (Mann-Whitney U-test or Wilcoxon signed rank test)
%   with significance level ALPHA (default is 0.01) and returns p values (WP).
%
%   DATA: Cell array or matrix with data for each group.
%   LABELS: Cell array with group labels.
%   ALPHA: Significance level (default is 0.01).
%   STR: 'nonpaired' for Mann-Whitney U-test, 'paired' for Wilcoxon signed rank test.
%   T: Title of the plot.
%   C: Cell array of colors for each group.
%
%   Example:
%       data = {randn(20, 1), randn(20, 1)};
%       labels = {'Group A', 'Group B'};
%       [H, Wp] = meanSE_barstat3(data, labels, 0.05, 'nonpaired', 'Comparison of Groups', {'r', 'b'});
%
%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025

% Input argument check
narginchk(2, 6);
if nargin < 4
    str = 'nonpaired';
end
if nargin < 3 || isempty(alpha)
    alpha = 0.01;   % default significance level
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
    bar(position, means(i), 'FaceColor', 'none', 'EdgeColor', c{i}, 'LineWidth', 1);
    errorbar(position, means(i), SEs(i), 'k', 'LineStyle', 'none', 'CapSize', 0);
    barPositions(i) = position; % Store the x-position for the bar
end
title(t); % Plot title

% Perform pairwise statistical tests
Wp = zeros(numGroups);
Wh = false(numGroups);
for i = 1:numGroups
    for j = i+1:numGroups
        switch str
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
