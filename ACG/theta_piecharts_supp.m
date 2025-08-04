function theta_piecharts_supp(ResponseCategoriCtrl, ResponseCategoriExp, TIC, TIE, resdir)
%THETA_PIECHARTS_SUPP Plots theta-group-based pie charts for Figure S12.
%
%   theta_piecharts_supp(ResponseCategoriCtrl, ResponseCategoriExp, TIC, TIE, resdir)
%   plots the proportion of strongly, moderately, and non-theta rhythmic
%   neurons across response categories for both control and experimental
%   datasets.
%
%   INPUTS:
%       ResponseCategoriCtrl - Cell array of response categories (control)
%       ResponseCategoriExp  - Cell array of response categories (experiment)
%       TIC                  - Theta index values for control neurons
%       TIE                  - Theta index values for experimental neurons
%       resdir               - Directory path to save output figures
%
%   See also: ACGMOD, CHISQUARETEST
%
%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025

    % Define group labels and plot appearance
    group_names = {'Inh','Act','Inh-Act','Act-Inh','NonResp'};
    colors_groups = {"#66CDAA","#EA9782","#0072BD", "#FF8383","#6161a7"};
    colors_theta = {"#EBDB34","#3495EB","#A07BE0"};
    explode5 = ones(1,5);
    explode3 = ones(1,3);
    theta_threshold_high = 0.5;
    theta_threshold_low = 0.1;

    % Process Experimental and Control datasets
    [~, thetaCellsExp, thetaProportionsExp] = processThetaGroups( ...
        ResponseCategoriExp, TIE, theta_threshold_high, theta_threshold_low, group_names);
    [~, thetaCellsCtrl, thetaProportionsCtrl] = processThetaGroups( ...
        ResponseCategoriCtrl, TIC, theta_threshold_high, theta_threshold_low, group_names);

    % Run chi-square test on group-level theta counts
    [Chi, pValue, dF] = chiSquareTest([thetaCellsExp; thetaCellsCtrl], 0.05); %#ok<ASGLU>

    % Plot the grouped pie charts
    plotThetaPieCharts(thetaCellsExp, thetaCellsCtrl, ...
        thetaProportionsExp, thetaProportionsCtrl, ...
        group_names, colors_groups, colors_theta, explode5, explode3);

    % Save figure to disk
    set(gcf, 'Renderer', 'painters');
    saveas(gcf, fullfile(resdir, 'Theta groups - S7.jpg'));
    saveas(gcf, fullfile(resdir, 'Theta groups - S7.svg'));
end

function [ratios, thetaCountsPerGroup, thetaProportions] = processThetaGroups(groups, thetaIndex, highTh, lowTh, group_names)
% Computes theta proportions within each response category.

    % Categorize each neuron based on theta index
    thetaGroup = categorizeTheta(thetaIndex, highTh, lowTh);

    % Logical index of all rhythmic neurons (strong + moderate)
    thetaRhythmic = thetaIndex >= lowTh;

    % Count total number of cells per group
    totalPerGroup = arrayfun(@(g) sum(strcmp(groups, g)), group_names);

    % Count number of theta-rhythmic neurons per group
    thetaCountsPerGroup = arrayfun(@(g) sum(strcmp(groups(thetaRhythmic), g)), group_names);

    % Compute theta ratio per group
    ratios = thetaCountsPerGroup ./ totalPerGroup;

    % Compute proportions of theta strength subgroups within each response category
    numThetaGroups = 3;
    thetaProportions = zeros(5, numThetaGroups);
    for g = 1:5
        idx = strcmp(groups, group_names{g});
        for t = 1:numThetaGroups
            thetaProportions(g, t) = sum(thetaGroup(idx) == t);
        end
        thetaProportions(g,:) = thetaProportions(g,:) / sum(thetaProportions(g,:));
    end
end

function thetaGroup = categorizeTheta(thetaIndex, highTh, lowTh)
% Categorizes neurons into theta strength groups.

    thetaGroup = zeros(size(thetaIndex));
    thetaGroup(thetaIndex > highTh) = 1;                               % Strongly theta-rhythmic
    thetaGroup(thetaIndex <= highTh & thetaIndex > lowTh) = 2;        % Moderately theta-rhythmic
    thetaGroup(thetaIndex <= lowTh) = 3;                              % Non-theta rhythmic
end

function plotThetaPieCharts(expCounts, ctrlCounts, propExp, propCtrl, group_names, colors5, colors3, explode1, explode2)
% Generates tiled pie chart layout for theta analysis.
    
    figure;
    t = tiledlayout(4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

     % Big pies for total theta rhythmic cells (Experimental & Control)
    nexttile(t, [2 1]);  % Experimental - spans 2 rows, 1 column
    plotColoredPie(expCounts, explode1, colors5);
    title('DR-2AFC WM');
    legend(group_names, "Box", "off", 'Position',[0.006250001064369,0.479412181270961,0.162499997871263,0.169047614435355]);

    % Small pies for each group (Experimental - top row)
    for g = 1:5
        if g <= 3 
            tile = g + 1;
        else
            tile = g + 2;
        end
        nexttile(t, tile);  % Tiles 2 to 6
        plotColoredPie(propExp(g,:), explode2, colors3);
        title({group_names{g}; ''});
    end

    legend({'Strongly Theta rhythmic','Moderately Theta','Non-Theta rhythmic'}, ...
                 'Box','off', 'Position',[0.718802088877806,0.583134923078241,0.273214280286006,0.104761901994546]);

    nexttile;
    axis off;

    nexttile(t, [2 1]);  % Control (positioned in 3rd row, 1st col)
    plotColoredPie(ctrlCounts, explode1, colors5);
    title('DR-2AFC Control');
    legend(group_names, "Box", "off", 'Position',[0.00547619047619,-0.003273809637342,0.162499997871263,0.169047614435355]);

    % Small pies for each group (Control - bottom row)
    for g = 1:5
        if g <= 3 
            tile = g + 9;
        else
            tile = g + 10;
        end
        nexttile(t, tile);  % Tiles 12 to 16
        plotColoredPie(propCtrl(g,:), explode2, colors3);
        title({group_names{g}; ''});
    end

    legend({'Strongly Theta rhythmic','Moderately Theta rhythmic','Non-Theta rhythmic'}, ...
                 'Box','off', 'Position',[0.718802088877806,0.045039684983,0.273214280286005,0.104761901994546]);

end

function plotColoredPie(data, explode, colors)
% Draw a pie chart with custom colors.

    h = pie(data, explode);
    for i = 1:length(data)
        h(2 * i - 1).FaceColor = colors{i};
    end
end

