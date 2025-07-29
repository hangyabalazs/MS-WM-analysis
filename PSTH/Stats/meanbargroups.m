function [ax, lgd, Wp, means, SEs] = meanbargroups(data, l, alpha, txt, parts)
% MEANBARGROUPS Create grouped bar plots with statistical significance.
%
%   [AX, LGD, WP, MEANS, SEs] = MEANBARGROUPS(DATA, L, TXT, PARTS) creates a
%   grouped bar plot of the data in DATA, with labels L, and group labels TXT.
%   It performs nonparametric Mann-Whitney U-tests between paired bars and
%   returns p-values, means, and standard errors.
%
%   Inputs:
%       DATA    - Cell array containing data vectors for each bar (10 elements).
%       L       - Cell array of x-axis labels for each bar (10 elements).
%       ALPHA   - Significance level  for the statistical test
%       TXT     - Title for the plot.
%       PARTS   - Flag to select color scheme:
%               0 - Red and Blue (e.g., WM vs Ctrl)
%               1 - Green and Red (e.g., Correct vs Incorrect)
%
%   Outputs:
%       AX      - Axes handle of the plot.
%       LGD     - Legend handle.
%       WP      - Vector of p-values from Mann-Whitney U-tests.
%       MEANS   - Vector of mean values for each bar.
%       SEs     - Vector of standard errors for each bar.
%
%       meanbargroups(data, labels, title_str, 0);
%
%   See also: MEANSE_BAR
%
%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025

    % Create the bar plot with error bars and statistical test
    
    % Define colors for the boxes
    if parts==1
        colors = {'g','r','g','r','g','r','g','r','g','r'};
    else
        colors = {'r','b','r','b','r','b','r','b','r','b'};
    end
    
    % Create labels
    group_labels = {'Inh', 'Act', 'Inh-Act', 'Act-Inh', 'NonResp'};
    spacing = 1.5;
    barPositions = zeros(10, 1);   % Store x-positions for bars
    for i = 1:10
        % Calculate position for each bar, alternating between no space for pairs and spacing between them
        groupIndex = ceil(i / 2);  % Group index (1, 2, 3, ...)
        barPositions(i) = (groupIndex - 1) * (2 + spacing) + (mod(i, 2) == 0); % Adjust positions for pairs
    end
    
    % Create the boxplot with customizations
    figure;
    [~, Wp, means, SEs] = meanSE_bar(data,l, alpha,'nonpaired',txt,colors);
    ylabel('Firing rate (Hz)');
    
    % Positioning the xticks: Place between pairs
    xtick_positions = zeros(5, 1);
    for i = 1:5
        % The x-position for each pair of bars is the average of the two bar positions
        xtick_positions(i) = mean([barPositions(2*i-1), barPositions(2*i)]);
    end
    
    % Define xtick labels (group labels)
    xtick_labels = group_labels;
    ax=gca;
    set(ax, 'XTick', xtick_positions, 'XTickLabel', xtick_labels);
    
    % Add custom legend using patches
    if parts==0
        p1 = patch(NaN, NaN, 'w', 'EdgeColor', 'red', 'LineWidth', 1);
        p2 = patch(NaN, NaN, 'w', 'EdgeColor', 'blue', 'LineWidth', 1);
        lgd=legend([p1, p2], {'WM', 'Ctrl'}, 'Location', 'best','Box','off');
    else 
        p1 = patch(NaN, NaN, 'w', 'EdgeColor', 'green', 'LineWidth', 1);
        p2 = patch(NaN, NaN, 'w', 'EdgeColor', 'red', 'LineWidth', 1);
        lgd=legend([p1, p2], {'Correct', 'Incorrect'}, 'Location', 'best','Box','off');
    end
end