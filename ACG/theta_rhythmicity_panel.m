function theta_rhythmicity_panel(resdir, cleaned_data, tranf_type)
%THETA_RHYTHMICITY_PANEL   Creates the panels with theta rhythmicity data.
%   This function generates Figure 5, showing the theta rhythmicity analysis
%   of experimental and control groups, including raster plots, PSTHs,
%   heatmaps, average response plots, pie charts, and statistical summaries.
%
%   Input parameters:
%       resdir      - Path to the directory for saving the figures and output.
%       cleaned_data - Struct containing cleaned cell data and associated metrics.
%       tranf_type   - Type of transformation applied to FR ('submean' or other).
%
%   Output:
%       Saves figures (in .svg and .jpg) and statistical results (in Excel).
%
%   See also ACGMOD, CHISQUARETEST, THETA_SUPP_PANELS,
%   THETA_PIECHARTS_SUPP.
%
%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025
    close all  % Close any open figures to start fresh

    % Extract relevant data from the cleaned_data structure
    cellidsC = cleaned_data.cellids.cellids_ctrl_cl;
    cellidsE = cleaned_data.cellids.cellids_exp_cl;
    normPSTH_exp = cleaned_data.spsth_data.normPSTH_exp;
    normPSTH_ctrl = cleaned_data.spsth_data.normPSTH_ctrl;
    normPSTH_exp_cl = cleaned_data.spsth_data.normPSTH_exp_cl;
    normPSTH_ctrl_cl = cleaned_data.spsth_data.normPSTH_ctrl_cl;
    spsthExp_raw = cleaned_data.spsth_data.spsthExp_raw;
    spsthCtrl_raw = cleaned_data.spsth_data.spsthCtrl_raw;
    time = cleaned_data.time;
    
    % Remove NaN rows from normPSTH and keep corresponding spsth
    spsthExp = spsthExp_raw(sum(isnan(normPSTH_exp),2) == 0, :);
    spsthCtrl = spsthCtrl_raw(sum(isnan(normPSTH_ctrl),2) == 0, :);
    
    % Load Theta Index data
    TIC = cleaned_data.ThetaIndex.ThetaIndexCtrl;
    TIE = cleaned_data.ThetaIndex.ThetaIndexExp;
    
    % Choose cellbase (ensure cellbase environment is set)
    choosecb('MS_WM_EXP_cellbase');
    
    % Define color scheme
    colors3 = {"#EBDB34", "#3495EB", "#A07BE0"};
    
    % Define example cell ids for example ACGs
    example_cells = {'NWM15_200312a_5.2', 'NWM5_190224a_3.1', 'NWM15_200226a_5.3'};
    
    % Define temporal parameters for delay and halves
    delay_start = 0.2; 
    firsthalf_end = 0.6; 
    delay_end = 1; 
    
    % Compute indices for different temporal windows
    delay_inx = find(time >= delay_start & time <= delay_end);
    firsthalf_inx = find(time >= delay_start & time <= firsthalf_end);
    secondhalf_inx = find(time >= firsthalf_end & time <= delay_end);
    
    % Define thresholds and labels for rhythmicity grouping
    highThreshold = 0.5;
    lowThreshold = 0.1;
    labels = {'Strongly theta-rhythmic cells', 'Moderately theta-rhythmic cells', 'Non-theta-rhythmic cells'};
    
    % Load speaker image used for plots
    speaker_image = imread([resdir '\Speaker_picture.png']);
    
    % Group control and experimental cells by theta index
    [highThetaIndicesC, countsC, sorted_psthsC, psthsCtrl, groupboundariesCtrl, t, avgCT, SEC, cC] = ...
        group_cells(TIC, highThreshold, lowThreshold, normPSTH_ctrl_cl, time, [delay_start, delay_end], cellidsC);
    [highThetaIndicesE, countsE, sorted_psthsE, psthsExp, groupboundariesExp, t, avgET, SEE, cE] = ...
        group_cells(TIE, highThreshold, lowThreshold, normPSTH_exp_cl, time, [delay_start, delay_end], cellidsE);
    
    ThetaInx = zeros(1,3);
    acg_tpos = {200, 40, 34};  % Positioning for theta index text
    
    % Top panel: autocorrelograms and theta index annotation
    for iC = 1:3
        [~, ~, ThetaInx(1,iC)] = acgmod(example_cells{iC}, 0.8, 'isvisual', true, 'dt', 0.001);  % ACG & compute theta index
        x_lim = xlim;
        y_lim = ylim;
        ylabel('Count');
        hold on;
        % Display theta index value on the plot
        text(x_lim(2)-1500, y_lim(2)-acg_tpos{iC}, mat2str(round(ThetaInx(iC),3)), 'Color', 'black', 'FontSize', 8);
        set(gca, 'TickDir', 'out', 'Box', 'off');
        hold off;
    end
    
    % Initialize main figure
    figure('Units', 'centimeters', 'Position', [0, 0, 16, 20]);
    
    % Copy autocorrelogram subplots into main figure (top row)
    for i = 1:3
        figure(i)
        h = get(gcf, 'Children');
        newh = copyobj(h(1), 4);
        set(newh, 'Position', [0.05 + (i-1)*0.33, 0.8, 0.28, 0.165]);
    end
    
    % Experimental heatmap
    figure(4)
    subplot('Position', [0.07, 0.55, 0.3, 0.2]);
    plot_heatmap(t, psthsExp, groupboundariesExp, 30, speaker_image, 0.12, 0.735);
    
    % Experimental average PSTH
    subplot('Position', [0.43, 0.55, 0.2, 0.2]);
    plot_avg(t, avgET, SEE, colors3, speaker_image, 0.482, 0.735);
    
    % Experimental group distribution (pie chart)
    subplot('Position', [0.745, 0.52, 0.22, 0.22]);
    plot_piechart(countsE, 'Exp', colors3, labels);
    
    % Control heatmap
    subplot('Position', [0.07, 0.28, 0.3, 0.2]);
    plot_heatmap(t, psthsCtrl, groupboundariesCtrl, 30, speaker_image, 0.12, 0.465);
    
    % Control average PSTH
    subplot('Position', [0.43, 0.28, 0.2, 0.2]);
    plot_avg(t, avgCT, SEC, colors3, speaker_image, 0.482, 0.465);
    
    % Control group distribution (pie chart)
    subplot('Position', [0.745, 0.2, 0.22, 0.22]);
    plot_piechart(countsC, 'Control', colors3, labels);
    
    % Statistical test: Chi-square on group distributions
    [chi_stat, p_value, df] = chiSquareTest([countsE; countsC], 0.001);
    res = [chi_stat, p_value, df];
    
    % Bar plots and statistics on firing rate
    [M1, S1, M2, S2, Wp12, M3, S3, M4, S4, Wp34, M5, S5, M6, S6, Wp56] = ...
        plot_bars(tranf_type, normPSTH_exp_cl, highThetaIndicesE, delay_inx, ...
        normPSTH_ctrl_cl, highThetaIndicesC, spsthExp, spsthCtrl, ...
        firsthalf_inx, secondhalf_inx);
    
    % Compile and format statistical results
    if strcmp(tranf_type, 'submean')
        emptyC = cell(11,3);
        emptyC(:) = {''};
        results = [
            {'Chi-square Statistic', 'pvalue', 'Degrees of Freedom'};
            res;
            {'Theta groups', '', ''};
            string(labels);
            {'Exp', '', ''};
            countsE;
            {'Ctrl', '', ''};
            countsC;
            emptyC;
            {'Subtracted mean FR - SE - p', 'Stronlgy theta-rythmic cells', 'Stats for delay'};
            {M1, S1, ''};
            {M2, S2, Wp12};
            {'Subtracted mean FR - SE - p', 'Stronlgy theta-rythmic cells', 'Stats for 1st half delay'};
            {M3, S3, ''};
            {M4, S4, Wp34};
            {'Subtracted mean FR - SE - p', 'Stronlgy theta-rythmic cells', 'Stats for 2nd half delay'};
            {M5, S5, ''};
            {M6, S6, Wp56}
        ];
    else
        results = [
            {'Chi-square Statistic', 'pvalue', 'Degrees of Freedom'};
            res;
            {'Theta groups', '', ''};
            string(labels);
            {'Exp', '', ''};
            countsE;
            {'Ctrl', '', ''};
            countsC;
            {'', '', ''};
            {'Absolute mean FR - SE - p', 'Stronlgy theta-rythmic cells', 'Stats for delay'};
            {M1, S1, ''};
            {M2, S2, Wp12};
            {'Absolute mean FR - SE - p', 'Stronlgy theta-rythmic cells', 'Stats for 1st half delay'};
            {M3, S3, ''};
            {M4, S4, Wp34};
            {'Absolute mean FR - SE - p', 'Stronlgy theta-rythmic cells', 'Stats for 2nd half delay'};
            {M5, S5, ''};
            {M6, S6, Wp56};
            {'', '', ''};
        ];
    end
    
    % Save results to Excel
    filename = fullfile(resdir, 'StatResults.xlsx');
    writematrix(results, filename, 'Sheet', 4);
    
    % Copy bottom bar plots into figure
    for i = 5:7
        figure(i)
        h = get(gcf, 'Children');
        newh = copyobj(h(2), 4);
        set(newh, 'Position', [0.06 + (i-5)*0.33, 0.03, 0.28, 0.165]);
    end
    
    % Finalize and save figure
    figure(4);
    set(figure(4), 'Renderer', 'painters');
    fnm1 = [resdir '\' 'Fig7_Theta rhythmicity panel' tranf_type '.svg'];
    saveas(figure(4), fnm1)
    fnm2 = [resdir '\' 'Fig7_Theta rhythmicity panel' tranf_type '.jpg'];
    saveas(figure(4), fnm2)
    
    % S6A - Generate theta rhythmicity supplementary panel for Experimental group
    theta_supp_panels('Exp', resdir, TIE, normPSTH_exp_cl, time);  % Generate and save supplementary figure for experimental
    close all  % Close figures after generating
    
    % S6B - Generate theta rhythmicity supplementary panel for Control group
    theta_supp_panels('Ctrl', resdir, TIC, normPSTH_ctrl_cl, time);  % Generate and save supplementary figure for control
    close all  % Close all figures to avoid overlap in following steps
    
    % S7 - Generate supplementary pie charts comparing delay response categories
    RC = cleaned_data.delay_response.ResponseCategoriCtrl;  % Extract response category data for control
    RE = cleaned_data.delay_response.ResponseCategoriExp;   % Extract response category data for experimental
    
    theta_piecharts_supp(RC, RE, TIC, TIE, resdir);  % Generate and save supplementary pie charts comparing groups

end 


function [highThetaIndices, counts, sorted_psths, psths, groupboundaries, time_cut, averages, SE, c]...
    = group_cells(TI, highThreshold, lowThreshold, normPSTH, time, delay, cellids)
% Groups cells based on theta index and prepares PSTHs and statistics for each group.

    % Find all neuron indices
    allIndices = 1:length(TI);
    
    % Categorize neurons based on Theta Index thresholds
    highThetaIndices = allIndices(TI > highThreshold);
    middleThetaIndices = allIndices(TI < highThreshold & TI > lowThreshold);
    lowThetaIndices = allIndices(TI < lowThreshold);
    
    % Count neurons per group
    numHighTheta = numel(highThetaIndices);
    numMiddleTheta = numel(middleThetaIndices);
    numLowTheta = numel(lowThetaIndices);
    counts = [numHighTheta, numMiddleTheta, numLowTheta];
    
    % Assign group labels: 1 = high, 2 = middle, 3 = low
    thetaGroup = zeros(size(TI));
    thetaGroup(highThetaIndices) = 1;
    thetaGroup(middleThetaIndices) = 2;
    thetaGroup(lowThetaIndices) = 3;
    
    % Preallocate for each group
    numGroups = 3;
    psths = cell(1, numGroups);
    averages = cell(1, numGroups);
    SE = cell(1, numGroups);
    
    % Extract PSTHs, compute means and SEs per group
    for g = 1:numGroups
        groupIndices = find(thetaGroup == g);
        psths{g} = normPSTH(groupIndices, :);
        averages{g} = mean(psths{g}, 1);
        SE{g} = std(psths{g})/sqrt(size(psths{g}, 1)); % Standard error
        c{g} = cellids(groupIndices); % Cell IDs for each group
    end
    
    % Sort PSTHs within each group by peak response in delay window
    sorted_psths = cell(1, numGroups);
    for g = 1:numGroups
        [mx, ~] = max(psths{g}(:, time >= delay(1) & time <= delay(2)), [], 2);
        [~, srtinx] = sort(mx, 'descend');
        sorted_psths{g} = psths{g}(srtinx, :);
    end
    
    % Combine sorted PSTHs into a single matrix
    psths = nan(length(TI), size(time,2));
    psths(1:counts(1), :) = sorted_psths{1};
    psths(counts(1)+1:counts(1)+counts(2), :) = sorted_psths{2};
    psths(counts(1)+counts(2)+1:end, :) = sorted_psths{3};
    
    % Store boundaries between groups for later visualization
    groupboundaries = [counts(1), sum(counts(1:2))];
    
    % Define time cut limits for trimming
    [~, time_start_index] = min(abs(time - (-0.504)));
    [~, time_end_index] = min(abs(time - 1.505));
    
    % Cut time vector and align all data accordingly
    time_cut = time(:, time_start_index:time_end_index);
    psths = psths(:, time_start_index:time_end_index);
    
    % Trim group statistics to the selected time window
    for iG = 1:numGroups
        SE{iG} = SE{iG}(:, time_start_index:time_end_index);
        averages{iG} = averages{iG}(:, time_start_index:time_end_index);
    end

end


function plot_heatmap(time, psths, groupboundaries, width, speaker_image, speaker_x, speaker_y)
% Plots a heatmap of PSTHs sorted by peak activity, with group boundaries and a speaker image overlay.

    % Display the PSTH matrix as a heatmap image
    imagesc(time, 1:size(psths, 1), psths);
    
    % Label the y-axis and set its limits to include padding on top and bottom
    ylabel('Neuron #');
    ylim([-width; size(psths, 1) + 0.5]);
    y = ylim;
    
    % Add colorbar to show scale of firing rate
    colorbar;
    clim([-9.51373008123223 9.51373008123223]); % Set color limits manually
    c = colorbar;
    c.Label.String = 'Normalized FR'; % Label the colorbar
    c.Label.Rotation = -90;           % Rotate label vertically
    c.TickDirection = 'out';          % Ticks face outward
    
    % Label x-axis
    xlabel('Time from cue onset (s)');
    
    hold on;
    
    % Plot horizontal lines to mark group boundaries
    for iC = 1:2
        line([time(:,1), time(:,end)], ...
             [groupboundaries(iC), groupboundaries(iC)], ...
             'Color', [0.5020 0.5020 0.5020], 'LineWidth', 0.5);
    end
    
    % Overlay speaker image at the specified position
    plot_fading_rectangles([1, size(psths, 1)], speaker_image, -width, width, speaker_x, speaker_y)
    
    hold off;
end


function plot_avg(time, average, SE, colors, speaker_image,  speaker_x,  speaker_y)
% Plots the average PSTH with shaded standard error bars for each neuron group.

    % Determine the number of theta rhythmicity groups
    numGroups = size(SE,2);
    
    % Plot each group's average PSTH with error shading
    for iC = 1:numGroups
        errorshade(time, average{iC}, SE{iC}, ...
            'LineColor', colors{iC}, 'ShadeColor', 'black', 'FaceAlpha', 0.4); % shaded error plot
    
        ylim([-3, 3]);                               % y-axis limits for firing rate
        xlim([time(:,1), time(:,end)]);              % x-axis based on time
        y = ylim;                                    % store y-limits (if needed)
        xlabel('Time from cue onset (s)');           % x-axis label
        ylabel('Average SPSTH');                     % y-axis label
    end 
    
    % Title based on condition
    if speaker_y == 0.735
        title('DR-2AFC WM');                         % Experimental condition
    else 
        title('DR-2AFC Control');                    % Control condition
    end
    
    % Overlay speaker image graphic
    hold on
    plot_fading_rectangles([-3;2.75], speaker_image, 2.75, 0.25, speaker_x, speaker_y);
    hold off;

end


function plot_fading_rectangles(y_limits, speaker_image, y_pos, height, speaker_x, speaker_y)
% Draws a fading effect over a stimulus period and overlays a speaker image.

    fade_length = 50; % Number of rectangles for fading transparency effect
    transparency_gradient = linspace(0.4, 0.02, fade_length); % Transparency gradient (opaque to almost transparent)
    x_pos = linspace(1, 1.5, fade_length); % Horizontal position of fading rectangles
    
    % Define main rectangle start and end
    start_index = 0.2;
    end_index = 1;
    
    % Color for solid and fading rectangles
    solid_color=[0.811764705882353, 0.250980392156863, 0.513725490196078]; % RGB solid fill
    fade_color = [0.811764705882353, 0.250980392156863, 0.513725490196078, 0.5]; % Same RGB, starting alpha
    
    % Plot vertical cue lines
    line([0, 0], [y_limits(1), y_limits(2)], 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1);
    if speaker_y == 0.735
        % Cue marker line
        line([0.2, 0.2], [y_limits(1), y_limits(2)], 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1);
    end
    
    % Solid bar showing stimulus duration
    rectangle('Position', [start_index, y_pos, end_index - start_index, height], ...
        'FaceColor', solid_color, 'EdgeColor', solid_color);

    % Loop to create a fade-out effect by stacking transparent rectangles
    for iC = 1:fade_length
        fade_color(4) = transparency_gradient(iC); % Update alpha value
        rectangle('Position', [x_pos(iC), y_pos, 0.01, height], ...
            'FaceColor', fade_color, 'EdgeColor', fade_color);
    end
    
    % Format current axis
    ax = gca;
    ax.TickDir = 'out'; % Ticks pointing outward
    ax.Box = 'off';     % Remove box frame
    
    % Create overlay axes for the speaker image
    ax2 = axes('Position', [speaker_x, speaker_y, 0.02, 0.02]); % Custom position
    imshow(speaker_image); % Show image
    axis('off');            % Hide axis ticks/labels
    grid on;

end


function plot_piechart(groupCounts, dataset, colors, labels)
% Creates a pie chart to visualize the proportion of neurons in different response categories.

    numGroups = size(groupCounts, 2);       % Number of groups/categories
    e = ones(1, numGroups);                 % Explode all slices equally
    h = pie(groupCounts, e);               % Create pie chart
    
    % Assign color to each slice
    for iC = 1:numGroups
        h(2 * iC - 1).FaceColor = colors{iC};  % Set face color of each pie slice
    end
    
    % Add legend only for control dataset
    if strcmp(dataset, 'Control')
        legend(labels, 'Position', [0.7307, 0.4385, 0.2, 0.1], 'Box', 'off');
    end
    
    hold off;

end


function [M1, S1, M2, S2, Wp12, M3, S3, M4, S4, Wp34, M5, S5, M6, S6, Wp56] = ...
    plot_bars(tranf_type, normPSTH_exp_cl, highThetaIndicesE, delay_inx, ...
              normPSTH_ctrl_cl, highThetaIndicesC, spsthExp, spsthCtrl, ...
              firsthalf_inx, secondhalf_inx)
% Plots bar comparisons of firing rates between experimental and control high-theta cells across different delay periods.

    % Statistical significance threshold
    alpha = 0.01;
    
    % Full delay period
    if strcmp(tranf_type, 'submean')
        [M1, S1, M2, S2, Wp12] = barmeanstat( ...
            mean(normPSTH_exp_cl(highThetaIndicesE, delay_inx), 2), ...
            mean(normPSTH_ctrl_cl(highThetaIndicesC, delay_inx), 2), 'Exp', 'Ctrl', alpha, 'nonpaired', 0);
        ylim([-3 1]);
    else
        [M1, S1, M2, S2, Wp12] = barmeanstat( ...
            mean(spsthExp(highThetaIndicesE, delay_inx), 2), ...
            mean(spsthCtrl(highThetaIndicesC, delay_inx), 2), 'Exp', 'Ctrl', alpha, 'nonpaired', 0);
    end
    title('Strongly theta-rhythmic cells - Delay');
    ylabel('Mean FR (Hz)');
    
    % First half of delay
    if strcmp(tranf_type, 'submean')
        [M3, S3, M4, S4, Wp34] = barmeanstat( ...
            mean(normPSTH_exp_cl(highThetaIndicesE, firsthalf_inx), 2), ...
            mean(normPSTH_ctrl_cl(highThetaIndicesC, firsthalf_inx), 2), 'Exp', 'Ctrl', alpha, 'nonpaired', 0);
        ylim([-2.5 0.5]);
    else
        [M3, S3, M4, S4, Wp34] = barmeanstat( ...
            mean(spsthExp(highThetaIndicesE, firsthalf_inx), 2), ...
            mean(spsthCtrl(highThetaIndicesC, firsthalf_inx), 2), 'Exp', 'Ctrl', alpha, 'nonpaired', 0);
    end
    title('1st Half - Delay');
    
    % Second half of delay
    if strcmp(tranf_type, 'submean')
        [M5, S5, M6, S6, Wp56] = barmeanstat( ...
            mean(normPSTH_exp_cl(highThetaIndicesE, secondhalf_inx), 2), ...
            mean(normPSTH_ctrl_cl(highThetaIndicesC, secondhalf_inx), 2), 'Exp', 'Ctrl', alpha, 'nonpaired', 0);
        ylim([-3.5 1]);
    else
        [M5, S5, M6, S6, Wp56] = barmeanstat( ...
            mean(spsthExp(highThetaIndicesE, secondhalf_inx), 2), ...
            mean(spsthCtrl(highThetaIndicesC, secondhalf_inx), 2), 'Exp', 'Ctrl', alpha, 'nonpaired', 0);
    end
    title('2nd Half - Delay');

end
