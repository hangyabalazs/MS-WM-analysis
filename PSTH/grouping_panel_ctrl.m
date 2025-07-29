function grouping_panel_ctrl(resdir, speaker_dir, response_categori_ctrl, normPSTH_ctrl_cl)
%GROUPING_PANEL_CTRL Generates panel of plots for visualizing grouped neural responses.
%
%   GROUPING_PANEL_CTRL(RES_DIR, SPEAKER_DIR, RESPONSE_CATEGORIES, NORM_PSTH) generates a
%   multi-panel figure for Ctrl neurons containing:
%       - Heatmaps of normalized PSTHs for each group
%       - Grouped average PSTHs with error shading
%       - Raster and PSTH plots for example cells
%
%   Inputs:
%       RESDIR              - Directory to save output figures
%       SPEAKER_DIR         - Directory containing speaker image
%       RESPONSE_CATEGORI_CTRL - Cell array of response categories for each
%       neuron in Ctrl
%       NORMPSTH_CTRL_CL      - Matrix of normalized PSTH data (neurons x time)
%
%
%   See also: ULTIMATE_PSTH_WM.
%
%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025

    if isempty(response_categori_ctrl) || isempty(normPSTH_ctrl_cl)
        error('Input data is empty or invalid.');
    end
    
    close all
    
    % Choose ctrl cellbase
    choosecb('MS_WM_CTRL_cellbase');
    
    % Define example cells 
    excellids={'NWM19_200928b_1.2'	'NWM19_200924b_4.1'	'NWM19_200908a_2.1'	'NWM19_200914b_4.2'	'NWM19_200914a_1.2'};
    
    % Define parameters
    group_names={'Inh','Act','Inh-Act','Act-Inh', 'NonResp'};
    full_group_names={'Inhibited then activated','Activated then inhibited','Inhibited','Activated','Non-responsive'};
    numGroups=5;
    time=linspace(-2,6,8001);
    delay_start=0.2;
    delay_end=1;
    colors={"#66CDAA","#EA9782","#0072BD", "#FF8383","#6161a7"};
    
    % Raster + PSTH aligned to stimulus onset parameters
    alignevent = 'FixationBeginning';   % trigger event
    partition = 'all';   % partition trials
    wn = [-2 6];   % full raster window in seconds
    dt = 0.001;   % raster resolution, in seconds
    sigma = 0.02;   % controls smoothing for 'spsth'
    bwin = [-1.6 0];   % baseline window for MW-test
    twin = [0.2 1];   % test-window for MW-test
    
    % Process psth data
    % Group psths
    psths = cell(1, numGroups);
    for i = 1:length(group_names)
        groupLabel = group_names{i};
        % Find indices of cells belonging to this group
        groupIndices = find(strcmp(response_categori_ctrl, groupLabel));
        
        % Collect PSTHs for this group
        psths{i} = normPSTH_ctrl_cl(groupIndices,:);
    end
    
    % Sort psths
    sorted_psths = cell(1, numGroups);
    for g = 1:numGroups
        [mx, mxinx] = max(psths{g}(:, time >= delay_start & time <= delay_end), [], 2);
        [~, srtinx] = sort(mx, 'descend');
        sorted_psths{g} = psths{g}(srtinx, :);
    end
    
    % Compute average and SE for psths
    avg = cell(1,5);
    SE = cell(1,5);
    for g = 1:numGroups
        avg{g} = mean(psths{g}, 1);
        SE{g} = std(psths{g})/sqrt(size(psths{g}, 1));
    end
    
    % Build panel
    % 1st subplot 
    % Average
    figure(1)
    average_plot(time, avg{1}, SE{1}, colors{1});
    
    figure(2)
    scatter_plot(excellids{3}, alignevent, sigma, partition, wn);
    
    figure(3)
    PSTH_plot(excellids{3}, alignevent, wn, dt, sigma, partition, bwin, twin, [1, 0.3], [10, 2], [0, 40]);
    
    figure(5)
    subplot('Position', [0.075, 0.7673, 0.12, 0.165])
    heatmap(time, sorted_psths{1}, 25, {full_group_names{3};''}, 0.102, speaker_dir);
    
    % Assembly of plots in 1 figure
    figure(1)
    h = get(gcf, 'Children');
    newh = copyobj(h,5);
    set(newh, 'Position',[0.075,0.55,0.12,0.165])
    
    figure(2)
    h = get(gcf, 'Children');
    newh = copyobj(h(1),5);
    set(newh, 'Position', [0.075, 0.32, 0.12, 0.165]);
    
    figure(4)
    h = get(gcf, 'Children');
    newh = copyobj(h(1),5);
    set(newh, 'Position', [0.075, 0.12, 0.12, 0.165]);
    xlim([-0.5, 1.5]);
    close 1 2 3 4
    
    % 2nd subplot
    figure(1)
    average_plot(time, avg{2}, SE{2}, colors{2})
    
    figure(2)
    scatter_plot(excellids{4}, alignevent, sigma, partition, wn);
    
    figure(3)
    PSTH_plot(excellids{4}, alignevent, wn, dt, sigma, partition, bwin, twin, [0, 1], [0, 10],[0, 100]);
    % 
    figure(5)
    subplot('Position', [0.26, 0.7673, 0.12, 0.165])
    heatmap(time, sorted_psths{2}, 10, {full_group_names{4};''}, 0.287,speaker_dir);
    % 
    figure(1)
    h = get(gcf, 'Children');
    newh = copyobj(h,5);
    set(newh, 'Position',[0.26,0.55,0.12,0.165])
    
    figure(2)
    h = get(gcf, 'Children');
    newh = copyobj(h(1),5);
    set(newh, 'Position',[0.26,0.32, 0.12,0.165]);
    
    figure(4)
    h = get(gcf, 'Children');
    newh = copyobj(h(1),5);
    set(newh, 'Position',[0.26, 0.12, 0.12, 0.165]);
    xlim([-0.5, 1.5]);
    close 1 2 3 4
    
    % 3rd subplot
    figure(1)
    average_plot(time, avg{3}, SE{3}, colors{3})
    
    figure(2)
    scatter_plot(excellids{1}, alignevent, sigma, partition, wn);
    
    figure(3)
    PSTH_plot(excellids{1}, alignevent, wn, dt, sigma, partition, bwin, twin, [1, 1], [1, 4], [0, 20]);
    
    figure(5)
    subplot('Position',[0.445,0.7673,0.12,0.165])
    heatmap(time, sorted_psths{3}, 0.9, {'Inhibited then';'Activated'}, 0.472,speaker_dir);
    
    figure(1)
    h = get(gcf,'Children');
    newh = copyobj(h,5);
    set(newh, 'Position',[0.445,0.55,0.12,0.165])
    
    figure(2)
    h = get(gcf,'Children');
    newh = copyobj(h(1),5);
    set(newh,'Position',[0.445,0.32, 0.12,0.165]);
    
    figure(4)
    h = get(gcf,'Children');
    newh = copyobj(h(1),5);
    set(newh,'Position',[0.445, 0.12, 0.12, 0.165]);
    xlim([-0.5, 1.5]);
    close 1 2 3 4
    
    % 4th subplot 
    figure(1)
    average_plot(time, avg{4}, SE{4}, colors{4})
    
    figure(2)
    scatter_plot(excellids{2}, alignevent, sigma, partition, wn);
    
    figure(3)
    PSTH_plot(excellids{2}, alignevent, wn, dt, sigma, partition, bwin, twin,  [1, 1], [0.1, 0.8],[0, 5]);
    
    figure(5)
    subplot('Position',[0.63,0.7673,0.12,0.165])
    heatmap(time, sorted_psths{4}, 0.08, {'Activated then';'Inhibited'}, 0.657,speaker_dir);
    
    figure(1)
    h = get(gcf,'Children');
    newh = copyobj(h,5);
    set(newh, 'Position',[0.63,0.55,0.12,0.165])
    
    figure(2)
    h = get(gcf,'Children');
    newh = copyobj(h(1),5);
    set(newh,'Position',[0.63,0.32, 0.12,0.165]);
    
    figure(4)
    h = get(gcf,'Children');
    newh = copyobj(h(1),5);
    set(newh,'Position',[0.63, 0.12, 0.12, 0.165]);
    xlim([-0.5, 1.5]);
    close 1 2 3 4
    
    % 5th subplot
    
    figure(1)
    average_plot(time, avg{5}, SE{5}, colors{5})
    
    figure(2)
    scatter_plot(excellids{5}, alignevent, sigma, partition, wn);
    
    figure(3)
    PSTH_plot(excellids{5}, alignevent, wn, dt, sigma, partition, bwin, twin,  [0, 0], [0, 0], [0, 100]);
    
    figure(5)
    subplot('Position',[0.815,0.7673,0.2,0.165])
    heatmap(time, sorted_psths{5}, 10, {full_group_names{5};''}, 0.842,speaker_dir);
    
    figure(1)
    h = get(gcf,'Children');
    newh = copyobj(h,5);
    set(newh, 'Position',[0.815,0.55,0.12,0.165])
    
    figure(2)
    h = get(gcf,'Children');
    newh = copyobj(h(1),5);
    set(newh,'Position',[0.815,0.32, 0.12,0.165]);
    
    figure(4)
    h = get(gcf,'Children');
    newh = copyobj(h(1),5);
    set(newh,'Position',[0.815, 0.12, 0.12, 0.165]);
    xlim([-0.5, 1.5]);
    close 1 2 3 4
    
    % Saving
    set(figure(5), 'Renderer', 'painters');
    
    fnm = [resdir '\'  'FigS5C_ctrl.svg'];
    saveas(figure(5),fnm);
    fnmm = [resdir '\'  'FigS5C_ctrl.jpg'];
    saveas(figure(5),fnmm);
    
    close all
end 


function average_plot(time,avg, SE, clr)
% Plot group average PSTH

    errorshade(time(:,1451:3506),avg(:,1451:3506),SE(:,1451:3506), 'LineColor',clr , 'ShadeColor', 'black', 'FaceAlpha', 0.4);
    xlim([-0.5,1.5]);
    ylim( [-6,6]);
    ylabel('Average SPSTH');
    y=ylim;
    xticks([0,1]);
    ax=gca;
    ax.TickDir = 'out';
    line([0,0], [y(1),y(2)], 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1);
end

function scatter_plot(example_cell,alignevent,sigma,partition,wn)
% Plot raster plot for example cell

    viewcell2b(example_cell,'TriggerName',alignevent,'SortEvent',alignevent,'sigma',sigma,...
    'eventtype','behav','ShowEvents',{{alignevent}},'Partitions',partition,'window',wn,'PSTHPlot',false);
    handle = get(gcf,'Children');
    delete(handle(3));
    delete(handle(2));
    delete(handle(1));
    ax=handle(4);
    ax.XLim=[-0.5 1.5];
    ax.YAxisLocation='left';
    ax.YLabel.String='Trial #';
    ax.YLabel.VerticalAlignment='top';
    ax.TickDir = 'out';
    y=ylim;
    line([0,0], [y(1),y(2)], 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1);
end

function PSTH_plot(example_cell, alignevent, wn, dt, sigma, partition, bwin, twin,x_positions, y_positions, y_limits)
% Plot raster plot for example cell

    [~, stats, ~, ~] = ...
    ultimate_psth_wm(example_cell,'trial',alignevent,wn,...
    'dt',dt,'sigma',sigma,'parts',partition,'isadaptive',0,...
    'maxtrialno',Inf,'baselinewin',bwin,'testwin',twin,'relative_threshold',0.01,...
    'display',true, 'event_filter', 'lowfixation_wm'); % calculate psth
    xlim([-0.5, 1.5]);
    xlabel({'Time from';' cue onset (s)'});
    x=xlim;
    ylabel('Firing rate');
    ax2=gca;
    ax2.Box='off';
    ax2.TickDir = 'out';
    ylim(y_limits);
    xticks([0,1]);
    
    y=ylim;
    WPi=stats.Wpi;
    WPa=stats.Wpa;
    if ~isnan(WPi)
        try 
            formatted_WPi=format_p_value(WPi);
            text(x(1)+x_positions(1),y(2)-y_positions(1),formatted_WPi,'Color',[0.00,0.60,1.00],'HorizontalAlignment','center','FontSize',10);
        catch 
             text(x(1)+x_positions(1),y(2)-y_positions(1),mat2str(round(WPi,3)),'Color','white','HorizontalAlignment','center','FontSize',5);
        end 
    end
    
    if ~isnan(WPa)
        try 
            formatted_WPa=format_p_value(WPa);
            text(x(1)+x_positions(2),y(2)-y_positions(2),formatted_WPa,'Color',[1.00,0.00,0.00],'HorizontalAlignment','center','FontSize',10);
        catch 
            text(x(1)+x_positions(2),y(2)-y_positions(2),mat2str(round(WPa,3)),'Color','white','HorizontalAlignment','center','FontSize',5);
        end
    end
    
    line([0,0], y_limits, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1);

end

function heatmap(time,psths,y_limit1,tit,x_pos, speakerdir)
% Plot psth heatmap for group
    
    start_index = 0.2;
    end_index = 1;
    solid_color=[0.811764705882353, 0.250980392156863, 0.513725490196078]; %RGB color for delay bar
    fade_color = [0.811764705882353, 0.250980392156863, 0.513725490196078, 0.5]; % RGB color with transparency for fading bar
    fade_length = 50; % Length of fading effect
    transparency_gradient = linspace(0.4, 0.02, fade_length);
    x_position=linspace(1,1.5,fade_length);
    try
        speaker_image = imread(fullfile([speakerdir '\Speaker_picture.png']));
    catch ME
        disp('Error encountered:');
        disp(ME.message);
    end
    
    imagesc(time(:,1451:3506), 1:size(psths, 1), psths(:,1451:3506)); %heatmap
    clim([-9.51373008123223 9.51373008123223]); %same colorbar for all heatmaps
    xlim([-0.5 1.5]);
    if x_pos==0.657
        ylim([-y_limit1;size(psths,1)+0.6]);
    elseif x_pos==0.472
        ylim([-y_limit1;size(psths,1)+0.6]);
    else 
        ylim([-y_limit1;size(psths,1)]);
    end
    
    xticks([0,1]);
    ylabel('Neuron #');
    
    if x_pos==0.842
        colorbar;
        c = colorbar;
        c.Label.String = 'Normalized FR';
        c.Label.Rotation=-90;
        c.Label.VerticalAlignment='middle';
        c.TickDirection='out';
    end
    
    title(tit,'FontSize',8,'HorizontalAlignment','center');
    y=ylim;
    line([0,0], [1,y(2)], 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1);
    
    if x_pos==0.472
        for i = 1:fade_length
            fade_color(4) = transparency_gradient(i); % Update transparency for fading effect
            rectangle('Position', [x_position(i), -y_limit1, 0.01, y_limit1+0.55], 'FaceColor', fade_color,'EdgeColor',fade_color);
        end
        rectangle('Position', [start_index,-y_limit1,end_index-start_index,y_limit1+0.55], 'FaceColor', solid_color,'EdgeColor',solid_color);
    elseif x_pos==0.657
        for i = 1:fade_length
            fade_color(4) = transparency_gradient(i); % Update transparency for fading effect
            rectangle('Position', [x_position(i), -y_limit1-0.02, 0.01, y_limit1+0.482], 'FaceColor', fade_color,'EdgeColor',fade_color);
        end
        rectangle('Position', [start_index,-y_limit1-0.02,end_index-start_index,y_limit1+0.482], 'FaceColor', solid_color,'EdgeColor',solid_color);
    else 
        for i = 1:fade_length
            fade_color(4) = transparency_gradient(i); % Update transparency for fading effect
            rectangle('Position', [x_position(i), -y_limit1, 0.01, y_limit1], 'FaceColor', fade_color,'EdgeColor',fade_color);
        end
        rectangle('Position', [start_index,-y_limit1,end_index-start_index,y_limit1], 'FaceColor', solid_color,'EdgeColor',solid_color);
    end
    
    ax1=gca;
    ax1.TickDir = 'out';
    ax1.Box = 'off';
    ax2=axes('Position',[x_pos, 0.92, 0.018, 0.018]);
    imshow(speaker_image);
    axis('off', 'image');

end