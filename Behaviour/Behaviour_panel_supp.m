function Behaviour_panel_supp(xlsdir, resdir)
%BEHAVIOUR_PANEL_SUPP   Plots supplementary behavioral panels (Figs S1-2)
%   BEHAVIOUR_PANEL_SUPP(XLSDIR, RESDIR) generates supplementary behavioral
%   figures (Figures S1-2) showing performance trends, pre- and post-surgery
%   comparisons, and fixation time analysis for experimental and control mice.
%
%   Inputs:
%       xlsdir - Full path to the directory containing the Excel file with 
%                behavioral data.
%       resdir - Full path to the directory where output figures should be saved.
%
%   This function is  called from BEHAVIOUR_PANEL to generate the full
%   set of figures for behavioral data.
%
%   See also BEHAVIOUR_PANEL
%
%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025


% Close all open figures to avoid confusion
close all

% Define directories
datadir_pre_exp = [xlsdir '\Pre_Surgery_Performance_Exp'];
datadir_pre_ctrl = [xlsdir '\Pre_Surgery_Performance_Ctrl'];
choosecb('MS_WM_EXP_cellbase');
datadir_post_exp = getpref('cellbase','datapath');
choosecb('MS_WM_CTRL_cellbase');
datadir_post_ctrl = getpref('cellbase','datapath');

% Extract data from excel file
[mouse_pre_exp, mouse_post_exp, mouse_pre_ctrl, mouse_post_ctrl,...
   NWM19_fix, NWM21_fix, NWM22_fix] = extract_excel_data(xlsdir);

% Define Mouse IDs
mouse_ids_exp = {'NWM5','NWM7','NWM9','NWM15'};
mouse_ids_ctrl = {'NWM19', 'NWM21', 'NWM22'};

% Define colors
mouse_colors_correct_exp = linspace(1, 0.5, length(mouse_ids_exp))' * [0.4, 1, 0.4];
mouse_colors_incorrect_exp = linspace(1, 0.5, length(mouse_ids_exp))' * [1, 0.4, 0.4];

mouse_colors_correct_ctrl = linspace(1, 0.5, length(mouse_ids_ctrl))' * [0.4, 1, 0.4];
mouse_colors_incorrect_ctrl = linspace(1, 0.5, length(mouse_ids_ctrl))' * [1, 0.4, 0.4];

% FIGURE S1: Experimental Mice 
fig_w = 20;

figure('Units', 'centimeters', 'Position', [0, 0, fig_w, 12]);
tlo = tiledlayout(3, 7, 'TileSpacing', 'compact', 'Padding', 'compact');

% Store all line/bar handles for legend
allHandles = [];

% Row 1: Pre and Post Surgery Performance (Side by Side)
nexttile(tlo, [1, 3]);
hold on;
h1 = plotMouseData(mouse_pre_exp, mouse_colors_correct_exp, mouse_colors_incorrect_exp, mouse_ids_exp);
allHandles = [allHandles; h1(:)];
ylabel('Completed Trials (%)', 'FontSize', 10);
xlabel('Training Session');
hold off;

nexttile(tlo, [1, 3]);
hold on;
plotMouseData(mouse_post_exp, mouse_colors_correct_exp, mouse_colors_incorrect_exp, mouse_ids_exp);
xlabel('Training Session');
hold off;

nexttile;
axis off;

% Row 2 & 3: pre-vs-post analysis, probability density estimates + bar
% plots
allHandles = pre_vs_post_perf(datadir_pre_exp, datadir_post_exp, 'exp', mouse_ids_exp, allHandles);

% Legend Tile (Rightmost column)
lgd = legend(allHandles, 'Location', 'northwest', 'Box', 'off', 'FontSize', 8);
set(lgd, 'Position', [0.6, 0.3, 0.4, 0.6]);

% Final formatting
set(findall(gcf, 'Type', 'text'), 'FontSize', 9);
saveas(gcf, fullfile(resdir, 'S1exp.svg'));
saveas(gcf, fullfile(resdir, 'S1exp.jpg'));


% FIGURE S2: Control Mice 

figure('Units', 'centimeters', 'Position', [0, 0, fig_w, 16]);
tlo2 = tiledlayout(4, 7, 'TileSpacing', 'compact', 'Padding', 'compact');

allHandles_ctrl = [];

% Row 1: Pre and Post Surgery Performance
nexttile(tlo2, [1, 3]);
hold on;
h3 = plotMouseData(mouse_pre_ctrl, mouse_colors_correct_ctrl, mouse_colors_incorrect_ctrl, mouse_ids_ctrl);
allHandles_ctrl = [allHandles_ctrl; h3(:)];
ylabel('Completed Trials (%)', 'FontSize', 10);
xlabel('Training Session');
hold off;

nexttile(tlo2, [1, 3]);
hold on;
plotMouseData(mouse_post_ctrl, mouse_colors_correct_ctrl, mouse_colors_incorrect_ctrl, mouse_ids_ctrl);
xlabel('Training Session');
hold off;

% Row 2: Fixation Time Line & fixation time bar graph
nexttile(tlo2, [1, 3]); 
h5 = plot_fixationtime(NWM19_fix, NWM21_fix, NWM22_fix);
allHandles_ctrl = [allHandles_ctrl; h5(:)];

nexttile(tlo2, [1, 3]); 
h6 = plot_fixationbar(NWM19_fix, NWM21_fix, NWM22_fix);
allHandles_ctrl = [allHandles_ctrl; h6(:)];

% Row 3 & 4: pre-vs-post analysis, probability density estimates + bar
% plots
allHandles_ctrl = pre_vs_post_perf(datadir_pre_ctrl, datadir_post_ctrl, 'ctrl', mouse_ids_ctrl, allHandles_ctrl);

% --- Legend Tile (Rightmost column) ---
lgd = legend(allHandles_ctrl, 'Location', 'northwest', 'Box', 'off', 'FontSize', 8);
set(lgd, 'Position', [0.6, 0.3, 0.4, 0.6]);

% Final formatting
set(findall(gcf, 'Type', 'text'), 'FontSize', 9);
saveas(gcf, fullfile(resdir, 'S2ctrl.svg'));
saveas(gcf, fullfile(resdir, 'S2ctrl.jpg'));

end

function allHandles = pre_vs_post_perf(datadir_pre, datadir_post, txt, mice_names, allHandles)
    
    % Extract data
    [GvaluesBefore, BvaluesBefore, WvaluesBefore] = performance_extrapolation(datadir_pre, mice_names, 'pre');
    [GvaluesAfter, BvaluesAfter, WvaluesAfter] = performance_extrapolation(datadir_post, mice_names, 'post');
    
    % === Row 2: ksdensity plots ===
    nexttile([1 2]);
    h1 = before_vs_after_ksdensity({-WvaluesBefore,-WvaluesAfter,[GvaluesBefore BvaluesBefore], [GvaluesAfter BvaluesAfter]},...
        {[0.7, 0.7, 0.7],[0.5, 0.5, 0.5],[0.2, 0.8, 0.8],[0, 0.5, 1]}, 'a', txt);
    allHandles(end+1:end+length(h1)) = h1;
    
    nexttile([1 2]);
    h2 = before_vs_after_ksdensity({GvaluesBefore, GvaluesAfter}, {[0, 1, 0],[0, 0.5, 0]}, 'b', txt);
    allHandles(end+1:end+length(h2)) = h2;
    
    nexttile([1 2]);
    h3 = before_vs_after_ksdensity({BvaluesBefore, BvaluesAfter},{[1, 0, 0],[0.5, 0, 0]}, 'c', txt);
    allHandles(end+1:end+length(h3)) = h3;

    nexttile;
    axis off;
     
    % === Row 3: bar plots ===
    nexttile ([1 2]);
    colours = {[0.7 0.7 0.7],[0.5 0.5 0.5]};
    h4 = before_vs_after_bar(colours,'a', WvaluesBefore, GvaluesBefore, WvaluesAfter, GvaluesAfter, BvaluesBefore, BvaluesAfter, txt);
    allHandles(end+1:end+length(h4)) = h4;
    
    nexttile([1 2]);
    colours = {[0 1 0],[0 0.7 0]};
    h5 = before_vs_after_bar(colours, 'b', GvaluesBefore, BvaluesBefore, GvaluesAfter, BvaluesAfter, txt);
    allHandles(end+1:end+length(h5)) = h5;
    
    nexttile([1 2]);
    colours = {[1 0 0],[0.7 0 0]};
    h6 = before_vs_after_bar(colours, 'c', BvaluesBefore, GvaluesBefore, BvaluesAfter, GvaluesAfter, txt);
    allHandles(end+1:end+length(h6)) = h6;

end

function h = plotMouseData(mouse_data, colors_correct, colors_incorrect, mouse_ids)
% Plot performance curves (completed trials vs sessions)

    h = gobjects(2, length(mouse_data)); % Preallocate: 2 lines per mouse (correct, error)

    % Loop through each mouse
    for mouse = 1:length(mouse_data)
        % Training days for each mouse
        num_training_days = size(mouse_data{mouse}, 2);
        training_days = 1:num_training_days;

        if size(mouse_data{mouse}, 1) == 2  % If there are 2 rows
            % No summing needed, directly use data
            correct_trials = mouse_data{mouse}(1, :);
            incorrect_trials = mouse_data{mouse}(2, :);
        elseif size(mouse_data{mouse}, 1) == 4  % If there are 4 rows
            % Sum the first two rows for correct trials and last two for incorrect
            correct_trials = sum(mouse_data{mouse}(1:2, :), 1);
            incorrect_trials = sum(mouse_data{mouse}(3:4, :), 1);
        else
            error('Unsupported data format for mouse %d', mouse);
        end
 
        % Plot correct trials with DisplayName
        h(1, mouse) = plot(training_days, correct_trials, 'Color', colors_correct(mouse, :), ...
            'LineWidth', 2, 'DisplayName', sprintf('%s Correct trials', mouse_ids{mouse}));
        hold on;

        % Plot incorrect trials
        h(2, mouse) = plot(training_days, incorrect_trials, 'Color', colors_incorrect(mouse, :), ...
            'LineWidth', 2, 'DisplayName', sprintf('%s Error trials', mouse_ids{mouse}));

    end

    % Labels and formatting
    xlabel('Training Session');
    ylabel('Completed Trials (%)', 'FontSize', 10);
    ylim([0, 1]);
    yticks([0, 1]);
    box off;
    set(gca, 'TickDir', 'out', 'TickLength', [0.015, 0.015]);
end

function [Gvalues, Bvalues, Wvalues] = performance_extrapolation(datadir,mice_id,tag)
% Extracts behavioral metrics from Bpod data : all time delta between central port 
% disengage and mandatory fixation end for completed trials or time missing to complete
% the mandatory fixation in aborted ones

    Gvalues = [];
    Bvalues = [];
    Wvalues = [];
    PIL = {'PokeInLight_1','PokeInLight_2','PokeInLight_3','PokeInLight_4','PokeInLight_5','PokeInLight_6'};
    
    for xx = 1:length(mice_id) % mouse loop
        Mouse_dir = dir(fullfile(datadir,mice_id{xx})); 
        for n = 3:length(Mouse_dir) % sessions loop
            try
                final = Mouse_dir(n).name;
                if strcmp(tag,'pre')
                    load([datadir '\' mice_id{xx} '\' final]);
                else
                    load([datadir '\' mice_id{xx} '\' final '\' mice_id{xx}  final '.mat' ]); % load bpod files
                end
                ntrials = SessionData.nTrials;
                Spans = zeros (1,ntrials);
                for nnn = 1:ntrials
                    fn = fieldnames(SessionData.RawEvents.Trial{1,nnn}.States);
                    for jj = 1 : length(fn) % catch fixation blocks number completed
                        if  strcmp(fn{jj}(1:3), 'Pok')
                            piln = str2num(fn{jj}(end));
                            vls = SessionData.RawEvents.Trial{1,nnn}.States.(PIL{piln});
                            if ~isnan(vls)
                                pil_end = vls(end);
                            end
                        end
                    end
                    
                    if ~isnan(SessionData.RawEvents.Trial{1,nnn}.States.ITI(1)) || ~isnan(SessionData.RawEvents.Trial{1,nnn}.States.pITI(1))% Completed Trials
                        if ~isnan(SessionData.RawEvents.Trial{1,nnn}.States.ITI(1)) 
                            ending = SessionData.RawEvents.Trial{1,nnn}.States.ITI(1); % successful trial end ts
                        else
                            ending = SessionData.RawEvents.Trial{1,nnn}.States.pITI(1); % unsuccessful trial end ts
                        end
                        if numel(SessionData.RawEvents.Trial{1,nnn}.Events.Port3Out) == 1 % proper port disengage after mandatory fixation ts
                            disengage = SessionData.RawEvents.Trial{1,nnn}.Events.Port3Out(1);
                        else
                            val = SessionData.RawEvents.Trial{1,nnn}.Events.Port3Out;
                            disengage = val(val > pil_end & val < ending);
                            disengage = disengage(1);
                        end
                        Spans(nnn) =  disengage - pil_end; % delta between light out and disengagement
                        if piln == 2 % adjustment to protocol dpendent fixation block number
                            ManFix(nnn) = (SessionData.RawEvents.Trial{1,nnn}.States.PokeInLight_2(end) - SessionData.RawEvents.Trial{1,nnn}.States.PokeInLight_1(1)); %  mandatory fixation total length
                        else
                            ManFix(nnn) = 3*(SessionData.RawEvents.Trial{1,nnn}.States.PokeInLight_2(end) - SessionData.RawEvents.Trial{1,nnn}.States.PokeInLight_1(1)); %  mandatory fixation total length
                        end
                        if ~isnan(SessionData.RawEvents.Trial{1,nnn}.States.ITI(1)) % Outcome coding (succ. or not)
                            Outcomes(nnn) = 1;
                        else
                            Outcomes(nnn) = 2;
                        end
                    elseif ~isnan(SessionData.RawEvents.Trial{1,nnn}.States.aITI(1)) || ~isnan(SessionData.RawEvents.Trial{1,nnn}.States.aaITI(1)) % Abort / withdrawl
                        Spans(nnn) = SessionData.RawEvents.Trial{1,nnn}.States.aITI(1)- SessionData.RawEvents.Trial{1,nnn}.States.PokeInLight_1(1);
                        Outcomes(nnn) = 3;
                        ManFix(nnn) = 0;
                    else
                        Outcomes(nnn) = 4;
                        ManFix(nnn) = 0;
                    end
                end
                
                ManFix = round(ManFix,2); % calculate fixation that should have been in withdrawed trials
                steps = unique(ManFix);
                if length(steps) > 2 % elongation
                    boundaries = zeros(1, length(steps));
                    for yy = 2 : length(steps)
                        [~, b] = find(ManFix == steps(yy));
                        boundaries(yy) =  b(end);
                    end
                    boundaries(end) = length(ManFix); % for Aborted last trials
                    ManFix2 = [];
                    for yy = 2 : length(boundaries)
                        block = ones(1,(boundaries(yy)-boundaries(yy-1)))*steps(yy);
                        ManFix2 = horzcat(ManFix2, block);
                    end
                    ManFix = ManFix2;
                else
                    ManFix = ones(1,ntrials)*steps(2);
                end
                GoodP = find(Outcomes == 1);
                yGLP = Spans(GoodP);
                Gvalues = horzcat(Gvalues, yGLP); % time values for successfull trials
                
                Bad = find(Outcomes == 2);
                yBLP = Spans(Bad); %#ok<*FNDSB>
                Bvalues = horzcat(Bvalues, yBLP); % time values for unsuccessfull trials
                
                Abort = find(Outcomes == 3);
                yWEP = ManFix(Abort) - Spans(Abort);
                yWEP = yWEP(yWEP>0); % fixing of "Fake Negative" trials
                Wvalues = horzcat(Wvalues, yWEP); % time values for aborted trials
            catch
                disp(['Error in handling ' final])
            end
            clearvars  Outcomes Spans ManFix
        end
    end
end

function h = before_vs_after_ksdensity(values_struct, line_colour, str, txt)
% Plots kernel density estimates : ksdensity distribution plot of 2 events time distributions
% (values_struct) before and after surgery, represented with distinct colours (line_colour)

    hold on;
    h = gobjects(1, length(values_struct));

    % Define DisplayName labels based on panel type
    switch str
        case 'a'
            displayNames = { ...
                'Aborted trials before surgery', ...
                'Aborted trials after surgery', ...
                'Completed trials before surgery', ...
                'Completed trials after surgery' };
        case 'b'
            displayNames = { ...
                'Correct trials before surgery', ...
                'Correct trials after surgery' };
        case 'c'
            displayNames = { ...
                'Error trials before surgery', ...
                'Error trials after surgery' };
        otherwise
            % Fallback names
            displayNames = cell(1, length(values_struct));
            for i = 1:length(values_struct)
                if mod(i,2)==1
                    displayNames{i} = sprintf('Before surgery (trace %d)', i);
                else
                    displayNames{i} = sprintf('After surgery (trace %d)', i);
                end
            end
    end

    % Plot each ksdensity with DisplayName
    for i = 1:length(values_struct)
        [f1, x1] = ksdensity(values_struct{i});
        h(i) = plot(x1, f1, 'Color', line_colour{i}, 'LineWidth', 2, ...
            'DisplayName', displayNames{i});
    end

    % Labels
    xlabel('Time from central port disengage (s)', 'FontSize', 8);
    ylabel('Probability', 'FontSize', 8.5);

    % Axis limits based on data type
    if strcmp(str,'a')
        if strcmp(txt, 'exp')
            xlim([-2.5, 4]);
            xticks([-2.5, -2, 0, 2, 4]);
            xticklabels({'', '-2', '0', '2', '4'});
            yticks([0, 2, 4]);
            ylim([0, 4]);
        else
            xlim([-1.5, 3]);
            xticks([-1.5, -1, 0, 1, 2, 3]);
            xticklabels({'', '-1', '0', '1', '2', '3'});
            yticks([0, 5, 10, 15]);
            ylim([0, 15]);
        end
    elseif strcmp(str,'b') || strcmp(str,'c')
        if strcmp(txt, 'exp')
            xlim([-1, 4]);
            xticks([-1, 0, 2, 4]);
            xticklabels({'', '0', '2', '4'});
            yticks([0, 1, 2, 2.5]);
            yticklabels({'0', '1', '2', ''});
            ylim([0, 2.5]);
        else
            xlim([-1, 3]);
            xticks([-1, 0, 1, 2, 3]);
            yticks([0, 3, 5]);
            yticklabels({'0', '3', '5'});
            ylim([0, 5]);
        end
    end

    set(gca, 'TickDir', 'out');
    hold off;
end

function h = before_vs_after_bar(bar_colour, str, MainBefore, SideBefore, MainAfter, SideAfter, SideBefore2, SideAfter2, txt)
% Plots bar graphs with error bars to compare performance before/after surgery

    % Calculate percentages
    switch nargin
        case 7 % Completed trials only
            total_before = numel(MainBefore) + numel(SideBefore);
            total_after = numel(MainAfter) + numel(SideAfter);
            bef = numel(MainBefore) / total_before;
            aft = numel(MainAfter) / total_after;
        case 9 % All trials
            total_before = numel(MainBefore) + numel(SideBefore) + numel(SideBefore2);
            total_after = numel(MainAfter) + numel(SideAfter) + numel(SideAfter2);
            bef = numel(MainBefore) / total_before;
            aft = numel(MainAfter) / total_after;
    end

    b_bar = [bef*100, aft*100];

    hold on;

    % Define DisplayName based on bar type
    switch str
        case 'a'
            name_before = 'Aborted trials before surgery';
            name_after  = 'Aborted trials after surgery';
            ylabel_str = 'All Trials (%)';
        case 'b'
            name_before = 'Correct trials before surgery';
            name_after  = 'Correct trials after surgery';
            ylabel_str = 'Completed Trials (%)';
        case 'c'
            name_before = 'Error trials before surgery';
            name_after  = 'Error trials after surgery';
            ylabel_str = 'Completed Trials (%)';
        otherwise
            name_before = 'Before';
            name_after  = 'After';
            ylabel_str = 'Performance (%)';
    end

    % Plot bars with DisplayName
    hb1 = bar(0, b_bar(1), 1, 'FaceColor', [1, 1, 1], ...
        'EdgeColor', bar_colour{1}, 'LineWidth', 2, ...
        'DisplayName', name_before);

    hb2 = bar(3, b_bar(2), 1, 'FaceColor', [1, 1, 1], ...
        'EdgeColor', bar_colour{2}, 'LineWidth', 2, ...
        'DisplayName', name_after);

    % Axis formatting
    xticks([0, 3]);
    xticklabels({''});
    ylim([0, max(60, max(b_bar)*1.1)]); % Dynamic ylim with floor at 60
    if strcmp(str,'a')
        xticklabels({''});
        if strcmp(txt,'exp')
            yticks([0,30,60]);
            yticklabels({'0','','60'});
            ty_pos = - 8;
            ylim([0,60]);
        else 
            yticks([0,10,20]);
            yticklabels({'0','','20'});
            ylim([0,20]);
            ty_pos = -2.7;
        end
        text(0,ty_pos,{'Aborted trials'; 'before surgery'},'HorizontalAlignment','center','FontSize',8);
        text(3,ty_pos,{'Aborted trials'; 'after surgery'},'HorizontalAlignment','center','FontSize',8);
    elseif strcmp(str,'b')
        ylim([0,90]);
        xticklabels({''});
        text(0,-12,{'Correct trials'; 'before surgery'},'HorizontalAlignment','center','FontSize',8);
        text(3,-12,{'Correct trials'; 'after surgery'},'HorizontalAlignment','center','FontSize',8);
        yticks([0,45,90]);
        yticklabels({'0','','90'});
    else 
        ylim([0,90]);
        xticklabels({''});
        text(0,-12,{'Error trials'; 'before surgery'},'HorizontalAlignment','center','FontSize',8);
        text(3,-12,{'Error trials'; 'after surgery'},'HorizontalAlignment','center','FontSize',8);
        yticks([0,45,90]);
        yticklabels({'0','','90'});
    end
    ylabel(ylabel_str, 'FontSize', 9);
    set(gca, 'TickDir', 'out', 'XLim', [-1, 4]);

    hold off;

    % Return bar handles for legend
    h = [hb1, hb2];
end

function h = plot_fixationtime(NWM19, NWM21, NWM22)
    % Line plot of fixation time trends

     % Fixation increase for Ctrl
    % Colors for mouse
    color19 = [0.678, 0.847, 0.902];
    color21 = [0, 0, 1];
    color22 = [0, 0, 0.545];
    
    hold on;
    h(1) = plot(NWM19, 'color', color19, 'LineWidth', 2, 'DisplayName', 'NWM19');
    h(2) = plot(NWM21, 'color', color21, 'LineWidth', 2, 'DisplayName', 'NWM21');
    h(3) = plot(NWM22, 'color', color22, 'LineWidth', 2, 'DisplayName', 'NWM22');
    plot([0, length(NWM22)], [0, 0], 'k');
    xlim([0, length(NWM22)]);
    ylim([0, 1]);
    yticks([-2,0,1,2]);
    ylabel('Fixation length (s)', 'FontSize', 9);
    xlabel('Training Session');
    set(gca, 'TickDir', 'out', 'TitleHorizontalAlignment', 'left', 'TickLength', [0.02, 0.02]);
    hold off;
end

function h = plot_fixationbar(NWM19, NWM21, NWM22)
% Bar plot comparing first and last day fixation times

    % Colors for mice
    color19 = [0.678, 0.847, 0.902]; % Light blue
    color21 = [0, 0, 1];             % Blue
    color22 = [0, 0, 0.545];         % Dark blue

    hold on;

    % Grouped bar for mean fixation time
    barmeans = [mean([NWM19(1), NWM21(1), NWM22(1)]), ...
                mean([NWM19(end), NWM21(end), NWM22(end)])];
    hb = bar([0, 3], barmeans, 0.6, 'FaceColor', 'none', ...
        'EdgeColor', [0, 0, 0.2], 'LineWidth', 2, ...
        'DisplayName', 'Mean across mice');

    % Individual mouse lines (first to last)
    hl1 = plot([0, 3], [NWM19(1), NWM19(end)], 'Color', color19, ...
        'LineWidth', 2, 'DisplayName', 'NWM19');
    hl2 = plot([0, 3], [NWM21(1), NWM21(end)], 'Color', color21, ...
        'LineWidth', 2, 'DisplayName', 'NWM21');
    hl3 = plot([0, 3], [NWM22(1), NWM22(end)], 'Color', color22, ...
        'LineWidth', 2, 'DisplayName', 'NWM22');

    % Reference line at zero
    plot([0, 3], [0, 0], 'k', 'LineWidth', 1);

    % Labels and formatting
    ylabel('Fixation length (s)', 'FontSize', 9);
    xticks([0, 3]);
    xticklabels({});
    xlim([-1.2, 4.2]);
    ylim([0, 1.5]);
    yticks([0, 1, 2]);

    % Session labels under x-axis
    text(0, -0.4, {'First'; 'training session'}, 'HorizontalAlignment', 'center', 'FontSize', 9);
    text(3, -0.4, {'Last'; 'training session'}, 'HorizontalAlignment', 'center', 'FontSize', 9);

    % Tick direction
    set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);

    hold off;

    % Return all relevant graphic handles for legend
    h = [hb, hl1, hl2, hl3];
end

function [mouse_pre_exp, mouse_post_exp, mouse_pre_ctrl, mouse_post_ctrl,...
   NWM19_fix, NWM21_fix, NWM22_fix] = extract_excel_data(xlsdir)
% Reads data from Excel sheets
    
    % Pre & post surgery perf data for exp animals
    mouse_pre_exp=cell(1, size(4, 2) );
    mouse_pre_exp{1} = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'B85:AE88');
    mouse_pre_exp{2} = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'B96:AE99');
    mouse_pre_exp{3} = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'B109:BX112');
    mouse_pre_exp{4} = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'B122:AJ125');
    
    mouse_post_exp=cell(1, size(4, 2) );
    mouse_post_exp{1} = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'CP19:EB22');
    mouse_post_exp{2} = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'CP26:EL29');
    mouse_post_exp{3} = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'CP34:EF37');
    mouse_post_exp{4} = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'CP43:FX46');
    
    % Pre & post surgery perf data for ctrl animals
    mouse_pre_ctrl=cell(1, size(3, 2) );
    mouse_pre_ctrl{1} = xlsread([xlsdir '\MiceTime.xlsx'], 2, 'B8:P9');
    mouse_pre_ctrl{2} = xlsread([xlsdir '\MiceTime.xlsx'], 2, 'B14:O15');
    mouse_pre_ctrl{3} = xlsread([xlsdir '\MiceTime.xlsx'], 2, 'B21:Q22');
    
    mouse_post_ctrl=cell(1, size(3, 2) );
    mouse_post_ctrl{1} = xlsread([xlsdir '\MiceTime.xlsx'], 2, 'B31:AZ32');
    mouse_post_ctrl{2} = xlsread([xlsdir '\MiceTime.xlsx'], 2, 'B37:AZ38');
    mouse_post_ctrl{3} = xlsread([xlsdir '\MiceTime.xlsx'], 2, 'B43:BA44');
    
    % Fixation times - Ctrl mice 
    NWM19_fix = xlsread([xlsdir '\MiceTime.xlsx'], 4, 'B1:AB1');
    NWM21_fix = xlsread([xlsdir '\MiceTime.xlsx'], 4, 'B5:AD5');
    NWM22_fix = xlsread([xlsdir '\MiceTime.xlsx'], 4, 'B8:AE7');
end