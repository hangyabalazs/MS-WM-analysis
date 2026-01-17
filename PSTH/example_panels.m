function example_panels(resdir,  group_assignments, normPSTH, grouping_type, cohort)
%EXAMPLE_PANELS Generates panel of plots for visualizing grouped neural responses.
%
%   Modified to merge:
%       - "Inh" + "Inh-Act"  "Inhibited"
%       - "Act" + "Act-Inh"  "Activated"
%       - "NonResp" unchanged
%
%   Panel: 4 rows  3 columns:
%       Row 1: Heatmaps (merged groups)
%       Row 2: Average PSTHs (merged groups)
%       Row 3: Raster plots (example cells from original Inh, Act, NonResp)
%       Row 4: PSTH plots (same example cells)
%
%   Example cells unchanged: picked from original categories.

%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025

    close all
    
    if isempty(group_assignments) || isempty(normPSTH)
        error('Input data is empty or invalid.');
    end
    
    % Choose cellbase and define example cells based on grouping type
    % (clustering or stat based grouping)
     if strcmp(cohort, 'exp')
         choosecb('MS_WM_EXP_cellbase');
         if strcmp(grouping_type, 'stat')
            excellids = {'NWM15_200221a_4.1', 'NWM15_200220b_5.3', 'NWM15_200223a_1.3'};
         else
             excellids = {'NWM15_200219a_2.1', 'NWM15_200219b_7.4', 'NWM15_200223a_1.1'};
         end
     else 
         choosecb('MS_WM_CTRL_cellbase');
         if strcmp(grouping_type, 'stat')
            excellids = {'NWM19_200908a_2.1', 'NWM19_200914b_4.2', 'NWM19_200914a_1.2'};
         else
            excellids={'NWM19_200908a_2.1',	'NWM19_200908b_5.1',	'NWM19_200910a_3.1'};
         end
     end

    % Define parameters
     if strcmp(grouping_type, 'stat')
        original_group_names = {'Inh','Act','Inh-Act','Act-Inh','NonResp'};
        group_names = {'Inhibited', 'Activated', 'Non-responsive'};
     else 
        group_names={'Cluster 1','Cluster 2','Cluster 3'};
     end

     numGroups = 3;
     colors = {"#66CDAA", "#EA9782", "#6161a7"};  % Colors for 3 merged groups
    time = linspace(-2,6,8001);
    delay_start = 0.2;
    delay_end = 1;

    % Raster + PSTH aligned to stimulus onset parameters
    alignevent = 'FixationBeginning';
    partition = 'all';
    wn = [-2 6];
    dt = 0.001;
    sigma = 0.02;
    bwin = [-1.6 0];
    twin = [0.2 1];

    psths = cell(1, numGroups);
    if strcmp(grouping_type, 'stat')
        % GROUP PSTHS (original 5 groups)
        psths_original = cell(1, 5);
        for iC = 1:5
            groupLabel = original_group_names{iC};
            groupIndices = find(strcmp(group_assignments, groupLabel));
            psths_original{iC} = normPSTH(groupIndices, :);
        end

    % MERGE GROUPS for Heatmap & Average PSTH
        % Inhibited = Inh + Inh-Act
        psths{1} = [psths_original{1}; psths_original{3}];
        % Activated = Act + Act-Inh
        psths{2} = [psths_original{2}; psths_original{4}];
        % Non-responsive = unchanged
        psths{3} = psths_original{5};
    else
        for clu = 1:3
            groupIndices = group_assignments == clu;
            psths{clu} = normPSTH(groupIndices, :);
        end
    end

    % Sort each merged group by max during delay
    sorted_psths = cell(1, numGroups);
    avg = cell(1, numGroups);
    SE = cell(1, numGroups);
    for g = 1:numGroups
        delay_mask = time >= delay_start & time <= delay_end;
        [mx, ~] = max(psths{g}(:, delay_mask), [], 2);
        [~, srtinx] = sort(mx, 'descend');
        sorted_psths{g} = psths{g}(srtinx, :);
        avg{g} = mean(psths{g}, 1);
        SE{g} = std(psths{g}) / sqrt(size(psths{g}, 1));
    end

    % CREATE PANEL 
    figure('Units', 'centimeter', 'Position',[0,0, 14.8167, 16]);
    tiledlayout(4, 3, 'Padding', 'loose', 'TileSpacing', 'loose');

    % ROW 1: HEATMAPS
    for col = 1:numGroups
        nexttile;
        heatmap_plot(time, sorted_psths{col}, group_names{col},  col==3, cohort);
        if col ==1 
            ylabel('Neuron #');
       end
    end

    % ROW 2: AVERAGE PSTHs
    for col = 1:numGroups
        nexttile;
        average_plot(time, avg{col}, SE{col}, colors{col}, cohort);
        if col == 1
            ylabel('Avg. Norm. FR');
        end
    end

    % ROW 3: RASTER PLOTS 
    for col = 1:numGroups
        figure;
        scatter_plot(excellids{col}, alignevent, sigma, partition, wn, cohort);
        if col == 1
            ylabel('Trial #');
        end
    end

    % ROW 4: PSTH PLOTS 
    % Define custom y_limits and text positions per example cell
    if strcmp(grouping_type, 'stat') 
        if strcmp(cohort, 'exp')
            y_limits_cell = {[15,40], [0,12], [0,1.5]};  % Inh, Act, NonResp
            x_positions_cell = {[1.1; 0], [0; 1], [1; 1]};
            y_positions_cell = {[10; 0], [0; 2], [0.1; 0.1]};
        else
            y_limits_cell = {[0,40], [0,100], [0,100]};  % Inh, Act, NonResp  from your original Ctrl code
            x_positions_cell = {[1, 1], [0, 0], [0, 0]};  % From your Ctrl PSTH_plot calls
            y_positions_cell = {[1, 4], [0, 0], [0, 0]};
        end
    else
        if strcmp(cohort, 'exp')
            y_limits_cell = {[0,24], [0,26], [0,8]};  
            x_positions_cell = {[1.1; 0], [1.1; 0], [1.1; 0]};
            y_positions_cell = {[5; 0], [5; 2], [2; 0]};
        else 
            y_limits_cell = {[5,40], [0,10], [2,16]};  
            x_positions_cell = {[1.1; 0], [1.1; 0], [1.1; 0]};
            y_positions_cell = {[20; 0], [2; 0], [2; 0]};
        end
    end
    
    for col = 1:numGroups
        figure;
        PSTH_plot(excellids{col}, alignevent, wn, dt, sigma, partition, bwin, twin, ...
                  x_positions_cell{col}, y_positions_cell{col}, y_limits_cell{col}, cohort);
        if col == 1
            ylabel('Firing rate (Hz)');
        end
    end
 
        % Raster inhibited
        figure(2)
        h = get(gcf,'Children');
        newhIS = copyobj(h(1),1); % scatter
        figure(1);
        ax = nexttile(7);
        axis off;
        pos = ax.Position;
        set(newhIS, 'Position', pos);
        
        % PSTH Inhibited
        figure(6)
        h = get(gcf,'Children');
        newhIP = copyobj(h(1),1); % psth
        figure(1);
        ax = nexttile(10);
        axis off;
        pos = ax.Position;
        set(newhIP, 'Position', pos);
        xlim([-0.5, 1.5]);
    
        % Raster Activated
        figure(3)
        h = get(gcf,'Children');
        newhAS = copyobj(h(1),1); % scatter
        figure(1);
        ax = nexttile(8);
        axis off;
        pos = ax.Position;
        set(newhAS, 'Position', pos);
        
        % PSTH activated
        figure(8)
        h = get(gcf,'Children');
        newhAP = copyobj(h(1),1); % psth
        figure(1);
        ax = nexttile(11);
        axis off;
        pos = ax.Position;
        set(newhAP, 'Position', pos);
        xlim([-0.5, 1.5]);
    
        % Raster non responsive
        figure(4)
        h = get(gcf,'Children');
        newhNS = copyobj(h(1),1); % scatter
        figure(1);
        ax = nexttile(9);
        axis off;
        pos = ax.Position;
        set(newhNS, 'Position', pos);
        
        % PSTH nonresponsive
        figure(10)
        h = get(gcf,'Children');
        newhNP = copyobj(h(1),1); % psth
        figure(1);
        ax = nexttile(12);
        axis off;
        pos = ax.Position;
        set(newhNP, 'Position', pos);
        xlim([-0.5, 1.5]);
    
    % --- SAVE ---
    set(figure(1), 'Renderer', 'painters');
    if strcmp(grouping_type, 'stat') 
        if strcmp(cohort, 'exp')
            fname = 'Fig4C_exp';
        else
            fname = 'Fig3C_ctrl';
        end
    else 
        if strcmp(cohort, 'exp')
            fname = 'FigS4C_exp';
        else 
            fname = 'FigS5_ctrl';
        end
    end
    saveas(figure(1), fullfile(resdir, [fname '.svg']));
    saveas(figure(1), fullfile(resdir, [fname '.jpg']));
    close all;
end

% -------------------------------------------------------------------------
function average_plot(time, avg, SE, clr, chrt)
    start_idx = 1451;  % Corresponds to t=-0.5s
    end_idx = 3506;    % Corresponds to t=1.5s
    errorshade(time(start_idx:end_idx), avg(start_idx:end_idx), SE(start_idx:end_idx), ...
               'LineColor', clr, 'ShadeColor', 'black', 'FaceAlpha', 0.4);
    xlim([-0.5, 1.5]);
    ylim([-6, 6]);
    xticks([0, 1]);
    ax = gca;
    ax.TickDir = 'out';
    y = ylim;
    line([0,0], y, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1);
    if strcmp(chrt,'exp')
        line([0.2,0.2], y, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1);
    end
end

function scatter_plot(example_cell, alignevent, sigma, partition, wn, chrt)
    viewcell2b(example_cell, 'TriggerName', alignevent, 'SortEvent', alignevent, 'sigma', sigma, ...
               'eventtype', 'behav', 'ShowEvents', {{alignevent}}, 'Partitions', partition, ...
               'window', wn, 'PSTHPlot', false);
    handle = get(gcf,'Children');
    delete(handle(3));
    delete(handle(2));
    delete(handle(1));
    ax=handle(4);
    % delete(ax.Children(1:3));  % Keep only raster
    ax.XLim = [-0.5, 1.5];
    ax.YAxisLocation = 'left';
    ax.TickDir = 'out';
    y = ylim;
    line([0,0], y, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1);
    if strcmp(chrt,'exp')
        line([0.2,0.2], y, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1);
    end
end

function PSTH_plot(example_cell, alignevent, wn, dt, sigma, partition, bwin, twin, x_positions, y_positions, y_limits, chrt)
    [~, stats, ~, ~] = ultimate_psth_wm(example_cell, 'trial', alignevent, wn, ...
        'dt', dt, 'sigma', sigma, 'parts', partition, 'isadaptive', 0, ...
        'maxtrialno', Inf, 'baselinewin', bwin, 'testwin', twin, 'relative_threshold', 0.01, ...
        'display', true, 'event_filter', 'lowfixation_wm');

    xlim([-0.5, 1.5]);
    ylim(y_limits);
    xlabel({'Time from'; 'cue onset (s)'});
    xticks([0,1]);

    ax = gca;
    ax.Box = 'off';
    ax.TickDir = 'out';

    WPi = stats.Wpi;
    WPa = stats.Wpa;
    x = xlim;
    y = ylim;

    if ~isnan(WPi)
        try
            formatted_WPi = format_p_value(WPi);
            text(x(1) + x_positions(1), y(2) - y_positions(1), formatted_WPi, ...
                 'Color', [0.00, 0.60, 1.00], 'HorizontalAlignment', 'center', 'FontSize', 10);
        catch
            text(x(1) + x_positions(1), y(2) - y_positions(1), mat2str(round(WPi,3)), ...
                 'Color', 'white', 'HorizontalAlignment', 'center', 'FontSize', 5);
        end
    end

    if ~isnan(WPa)
        try
            formatted_WPa = format_p_value(WPa);
            text(x(1) + x_positions(2), y(2) - y_positions(2), formatted_WPa, ...
                 'Color', [1.00, 0.00, 0.00], 'HorizontalAlignment', 'center', 'FontSize', 10);
        catch
            text(x(1) + x_positions(2), y(2) - y_positions(2), mat2str(round(WPa,3)), ...
                 'Color', 'white', 'HorizontalAlignment', 'center', 'FontSize', 5);
        end
    end

    line([0,0], y_limits, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1);
    if strcmp(chrt,'exp')
        line([0.2,0.2], y_limits, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1);
    end
end

function heatmap_plot(time, psths, title_str,  show_colorbar, chrt)
    start_idx = 1451;  % t = -0.5s
    end_idx = 3506;    % t = 1.5s
    % delay_start = 0.2;
    % delay_end = 1;
    % 
    % solid_color = [0.8118, 0.2510, 0.5137];  % RGB for delay bar
    % fade_color = [0.8118, 0.2510, 0.5137, 0.5];
    % fade_length = 50;
    % transparency_gradient = linspace(0.4, 0.02, fade_length);
    % x_position = linspace(1, 1.5, fade_length);

    imagesc(time(start_idx:end_idx), 1:size(psths,1), psths(:,start_idx:end_idx));
    clim([-9.5137, 9.5137]);  % Fixed for consistency
    xlim([-0.5, 1.5]);
    ylim([1, size(psths,1)]);
    xticks([0,1]);
    title(title_str, 'FontSize', 9, 'FontWeight', 'bold');

    y = ylim;
    line([0,0], [1,y(2)], 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1);
    if strcmp(chrt,'exp')
        line([0.2,0.2], [1,y(2)], 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1);
    end
    % % Solid delay bar
    % rectangle('Position', [delay_start, 1, delay_end-delay_start, size(psths,1)], ...
    %           'FaceColor', solid_color, 'EdgeColor', solid_color);
    % 
    % % Fading bar
    % for i = 1:fade_length
    %     fade_color(4) = transparency_gradient(i);
    %     rectangle('Position', [x_position(i), 1, 0.01, size(psths,1)], ...
    %               'FaceColor', fade_color, 'EdgeColor', fade_color);
    % end

    ax = gca;
    ax.TickDir = 'out';
    ax.Box = 'off';

    if show_colorbar
        cb = colorbar;
        cb.Label.String = 'Norm. FR';
        cb.Label.Rotation = -90;
        cb.Label.VerticalAlignment = 'middle';
        cb.TickDirection = 'out';
    end

    
end