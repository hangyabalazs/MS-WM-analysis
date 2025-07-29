function theta_supp_panels(dataset, resdir, ThetaIndexMatrix, norm_spsth, time)
%THETA_SUPP_PANELS   Creates supplementary panels for theta example cells (Figs. S9 & S10).
%   THETA_SUPP_PANELS(DATASET, RESDIR, THETAINDEXMATRIX, NORM_SPSTH, TIME)
%   generates multi-panel figures showing example neurons sorted by theta
%   rhythmicity (strong, moderate, non-theta) using raster plots, PSTHs, and
%   heatmaps.
%
%   INPUTS:
%       dataset           - 'Ctrl' or 'Exp', specifying the dataset
%       resdir            - Directory to save output figures
%       ThetaIndexMatrix  - Vector of theta indices for all neurons
%       norm_spsth        - Matrix of normalized spike density (neurons x time)
%       time              - Time vector corresponding to norm_spsth
%
%   See also: ACGMOD, ERRORSHADE, VIEWCELL2B, ULTIMATE_PSTH_WM.
%
%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025

    close all  % Close all existing figures
    
    % Dataset-specific setup
    if strcmp(dataset,'Ctrl')
        choosecb('MS_WM_CTRL_cellbase');
        example_cells = {'NWM19_200914a_1.2', 'NWM19_200915a_8.1',  'NWM19_200910a_5.3'};
        x_positions = {[1; 0], [1; 1.5],[1; 0]};
        y_positions = {[10; 0], [2; 1.5], [1; 0]};
        y_limits = {[35,75], [2,16],[0,4]};
        heatmap_spacing = {2, 10, 35};
    elseif strcmp(dataset, 'Exp')
        choosecb('MS_WM_EXP_cellbase');
        example_cells = {'NWM15_200312a_5.2', 'NWM5_190224a_3.1', 'NWM15_200226a_5.3'};
        x_positions = {[1; 0], [1; 0], [0; 1.5]};
        y_positions = {[10; 0], [2; 0], [0; 2]};
        y_limits = {[10,60], [10,90], [0,8]};
        heatmap_spacing = {1.3, 10, 40};
    end
    
    % Parameters
    colors3 = {"#EBDB34", "#3495EB", "#A07BE0"};  % Colors for each theta group
    alignevent = 'FixationBeginning';
    sigma = 0.02;
    dt = 0.001;
    wn = [-2 6];              % Window 
    bwin = [-1.6 0];          % Baseline window
    twin = [0.2 1];           % Task window
    delay_start = 0.2;        % Start of task window
    delay_end = 1;            % End of task window
    avg_y_limits = [-3 3];
    avg_x_limits = [-0.5 1.5];
    
    % Categorize neurons by theta index
    highThreshold = 0.5;
    lowThreshold = 0.1;
    allIndices = 1:length(ThetaIndexMatrix);
    
    highThetaIndices = allIndices(ThetaIndexMatrix > highThreshold);
    middleThetaIndices = allIndices(ThetaIndexMatrix < highThreshold & ThetaIndexMatrix > lowThreshold);
    lowThetaIndices = allIndices(ThetaIndexMatrix < lowThreshold);
    
    thetaGroup = categorizeTheta(ThetaIndexMatrix, highThreshold, lowThreshold);
    
    % Split and sort PSTHs by theta group
    numGroups = 3;
    psths = cell(1, numGroups);
    sorted_psths = cell(1, numGroups);
    
    for g = 1:numGroups
        groupIndices = find(thetaGroup == g);
        psths{g} = norm_spsth(groupIndices, :);
        
        % Sort each PSTH by max activity in delay window
        [mx, ~] = max(psths{g}(:, time >= delay_start & time <= delay_end), [], 2);
        [~, srtinx] = sort(mx, 'descend');
        sorted_psths{g} = psths{g}(srtinx, :);
    end
    
    % Loop through each theta group (strong, moderate, weak)
    for i = 1:3
        % Get example cell and group indices
        switch i
            case 1
                cell_idx = highThetaIndices;
                label = {'Strongly Theta'; 'rhythmic'};
            case 2
                cell_idx = middleThetaIndices;
                label = {'Moderately Theta'; 'rhythmic'};
            case 3
                cell_idx = lowThetaIndices;
                label = {'Non-Theta'; 'rhythmic'};
        end
        
        % Average PSTH
        figure(1)
        averagePlot(time, norm_spsth(cell_idx, :), colors3{i}, dataset);
        xlim(avg_x_limits); ylim(avg_y_limits);
        line([0, 0], avg_y_limits, 'Color', [0.929 0.694 0.125], 'LineWidth', 1);
        if strcmp(dataset,'Exp')
            line([0.2, 0.2], avg_y_limits, 'Color', [0.929 0.694 0.125], 'LineWidth', 1);
        end
    
        % Raster plot
        figure(2)
        scatterPlot(example_cells{i}, alignevent, sigma, 'all', wn, dataset);
    
        % PSTH plot
        figure(4)
        PSTHPlot(example_cells{i}, alignevent, wn, dt, sigma, 'all', ...
            bwin, twin, x_positions{i}, y_positions{i}, y_limits{i}, dataset);
    
        % Heatmap 
        if i == 1
            figure('Units', 'centimeters', 'Position', [0, 0, 10, 12]);
        else 
            figure(5)
        end
        if i == 3
            subplot('Position', [0.08 + 0.31 * (i-1), 0.7673, 0.3, 0.165])
        else 
            subplot('Position', [0.08 + 0.31 * (i-1), 0.7673, 0.2, 0.165])
        end
        plot_heatmap(resdir, time, sorted_psths{i}, heatmap_spacing{i}, label, dataset, i);
    
        % Compose Panel
        figure(1)
        h = get(gcf, 'Children');
        newh = copyobj(h, 5);
        set(newh, 'Position', [0.08 + 0.31 * (i-1), 0.55, 0.2, 0.165], 'XTick', [0 1]);
    
        figure(2)
        h = get(gcf, 'Children');
        newh = copyobj(h(1), 5);
        set(newh, 'Position', [0.08 + 0.31 * (i-1), 0.32, 0.2, 0.165], 'XTick', [0 1]);
    
        figure(3)
        h = get(gcf, 'Children');
        newh = copyobj(h(1), 5);
        set(newh, 'Position', [0.08 + 0.31 * (i-1), 0.12, 0.2, 0.165], 'XTick', [0 1]);
        xlim([-0.5, 1.5]);
    
        % Close temporary figures
        close(1); close(2); close(3); close(4);
    end
    
    % Save Final Composite Figure
    figure(5)
    set(gcf, 'Renderer', 'painters');
    
    svg_name = fullfile(resdir, [dataset ' Panel Examples.svg']);
    jpg_name = fullfile(resdir, [dataset ' Panel Examples.jpg']);
    
    saveas(gcf, svg_name);
    saveas(gcf, jpg_name);

end


function thetaGroup = categorizeTheta(thetaIndex, highTh, lowTh)
% Categorizes neurons into theta rhythmicity groups based on thresholds.

    % Initialize group labels with zeros
    thetaGroup = zeros(size(thetaIndex));
    thetaGroup(thetaIndex > highTh) = 1;     % Group 1: Strongly theta rhythmic
    thetaGroup(thetaIndex <= highTh & thetaIndex > lowTh) = 2;     % Group 2: Moderately theta rhythmic
    thetaGroup(thetaIndex <= lowTh) = 3;     % Group 3: Non-theta
end

%---------------------------------------------------------------------------

function averagePlot(time, data, clr, str)
% Plots the average normalized PSTH with error shading.

    % Compute the mean across neurons (rows)
    mn = mean(data, 1);

    % Compute standard error of the mean across neurons
     SE = std(data)/sqrt(size(data, 1));

    % Plot mean and shaded error using custom function
    errorshade(time, mn, SE, ...
        'LineColor', clr, ...
        'ShadeColor', 'black', ...
        'FaceAlpha', 0.4);

    % Set axis properties
    xlim([-0.5, 1.5]);
    ylabel('Average SPSTH');
    ax = gca;
    ax.TickDir = 'out';

    % Add alignment line at time = 0
    y = ylim;  % Get y-limits to draw full-height lines
    line([0, 0], y, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1);

    % Optionally add second line for delay onset (0.2s) in 'Exp' dataset
    if strcmp(str, 'Exp')
        line([0.2, 0.2], y, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1);
    end

end


function [WPa, WPi] = PSTHPlot(example_cell, alignevent, wn, dt, sigma, ...
                                partition, bwin, twin, x_positions, y_positions, y_limits, str)
% Plot PSTH for a single example cell and annotate Wilcoxon p-values.

    % Compute PSTH using specified parameters and extract stats
    [~, stats, ~, ~] = ultimate_psth_wm( ...
        example_cell, 'trial', alignevent, wn, ...
        'dt', dt, 'sigma', sigma, 'parts', partition, ...
        'isadaptive', 0, 'maxtrialno', Inf, ...
        'baselinewin', bwin, 'testwin', twin, ...
        'relative_threshold', 0.01, ...
        'display', true, ...
        'event_filter', 'lowfixation_wm');

    % Set axis limits and labels
    xlim([-0.5, 1.5]);
    xlabel({'Time from', 'cue onset (s)'}, 'FontSize', 5);
    ylabel('Firing rate');
    ax = gca;
    ax.Box = 'off';
    ax.TickDir = 'out';
    ylim(y_limits);
    x = xlim;
    y = ylim;

    % Extract and display p-value for inhibition (WPi)
    WPi = stats.Wpi;
    if ~isnan(WPi)
        try
            formatted_WPi = format_p_value(WPi);
            text(x(1) + x_positions(1), y(2) - y_positions(1), formatted_WPi, ...
                 'Color', [0.00, 0.60, 1.00], 'HorizontalAlignment', 'center', 'FontSize', 10);
        catch
            text(x(1) + x_positions(1), y(2) - y_positions(1), ...
                 mat2str(round(WPi, 3)), 'Color', 'white', ...
                 'HorizontalAlignment', 'center', 'FontSize', 5);
        end
    end

    % Extract and display p-value for activation (WPa)
    WPa = stats.Wpa;
    if ~isnan(WPa)
        try
            formatted_WPa = format_p_value(WPa);
            text(x(1) + x_positions(2), y(2) - y_positions(2), formatted_WPa, ...
                 'Color', [1.00, 0.00, 0.00], 'HorizontalAlignment', 'center', 'FontSize', 10);
        catch
            text(x(1) + x_positions(2), y(2) - y_positions(2), ...
                 mat2str(round(WPa, 3)), 'Color', 'white', ...
                 'HorizontalAlignment', 'center', 'FontSize', 5);
        end
    end

    % Draw vertical line at time = 0
    line([0, 0], y_limits, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1);

    % Optional vertical line at 0.2s for experimental dataset
    if strcmp(str, 'Exp')
        line([0.2, 0.2], y_limits, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1);
    end

end


function scatterPlot(example_cell, alignevent, sigma, partition, wn, str, n)
% Visualize spike raster for a single example cell.

    if nargin < 7
        n = 'all';
    end

    % Generate raster plot using viewcell2b
    viewcell2b(example_cell, ...
        'TriggerName', alignevent, ...
        'SortEvent', alignevent, ...
        'sigma', sigma, ...
        'eventtype', 'behav', ...
        'ShowEvents', {{alignevent}}, ...
        'ShowEventsColors', {{[0.9290 0.6940 0.1250]}}, ...
        'Partitions', partition, ...
        'window', wn, ...
        'PSTHPlot', false, ...
        'Num2Plot', n);

    % Clean up extra subplots
    h = get(gcf, 'Children');
    delete(h(3)); delete(h(2)); delete(h(1));
    
    % Format raster plot axis
    ax = h(4);
    ax.XLim = [-0.5, 1.5];
    ax.YAxisLocation = 'left';
    ax.YLabel.String = 'Trial #';
    ax.YLabel.VerticalAlignment = 'top';
    ax.TickDir = 'out';
    y = ylim;

    % Add vertical line for cue onset at 0.2s (experimental)
    if strcmp(str, 'Exp')
        line([0.2, 0.2], [y(1), y(2)], 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1);
    end
end


function plot_heatmap(resdir, time, psths, spacing, group, dataset, plot_num)
% Plot heatmap of normalized firing rates with annotation.

    % Delay bar timing
    start_index = 0.2;
    end_index = 1;

    % Bar colors
    solid_color=[0.811764705882353, 0.250980392156863, 0.513725490196078]; %RGB color for delay bar
    fade_color = [0.811764705882353, 0.250980392156863, 0.513725490196078, 0.5]; % RGB color with transparency for fading bar
    fade_length = 50;
    transparency_gradient = linspace(0.4, 0.02, fade_length);
    x_position = linspace(1, 1.5, fade_length);

    % Load speaker image
    speaker_image = imread([resdir '\Speaker_picture.png']);

    % Main heatmap plot
    imagesc(time, 1:size(psths, 1), psths);
    clim([-9.5137, 9.5137]); % Color range fixed for normalization
    xlim([-0.5, 1.5]);
    xticks([0 1]);
    ylabel('Neuron #');

    % Axis limits based on plot type
    if plot_num == 3
        colorbar;
        c = colorbar;
        c.Label.String = 'Normalized FR';
        c.Label.Rotation = -90;
        c.Label.VerticalAlignment = 'middle';
        c.TickDirection = 'out';
        ylim([-spacing, size(psths,1)]);
    elseif plot_num == 1
        ylim([-spacing, size(psths,1) + 1]);
    else
        ylim([-spacing, size(psths,1)]);
    end

    % Title and vertical event lines
    title(group, 'FontSize', 8, 'HorizontalAlignment', 'center');
    y = ylim;
    line([0, 0], [1, y(2)], 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1);

    if strcmp(dataset, 'Exp')
        line([0.2, 0.2], [1, y(2)], 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1);
    end

    % Solid rectangle for delay bar
    rectangle('Position', [start_index, -spacing, end_index - start_index, spacing], ...
              'FaceColor', solid_color, 'EdgeColor', solid_color);

    % Fade-out bar overlay
    for i = 1:fade_length
        fade_color(4) = transparency_gradient(i); % Update alpha
        rectangle('Position', [x_position(i), -spacing, 0.01, spacing], ...
                  'FaceColor', fade_color, 'EdgeColor', fade_color);
    end

    % Format axis
    ax1 = gca;
    ax1.TickDir = 'out';
    ax1.Box = 'off';

    % Place speaker icon
    if plot_num == 1
        speaker_x = 0.125;
    elseif plot_num == 2
        speaker_x = 0.435;
    else
        speaker_x = 0.735;
    end

    ax2 = axes('Position', [speaker_x, 0.915, 0.023, 0.023]);
    imshow(speaker_image);
    axis off;
end