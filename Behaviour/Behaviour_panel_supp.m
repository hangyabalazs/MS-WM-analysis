function Behaviour_panel_supp(xlsdir, resdir)
%BEHAVIOUR_PANEL_SUPP   Plots supplementary behavioral panels (Figs S1�S4)
%   BEHAVIOUR_PANEL_SUPP(XLSDIR, RESDIR) generates supplementary behavioral
%   figures (Figures S1-S4) showing performance trends, pre- and post-surgery
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

% Define directories containing pre- and post-surgery performance data for
% Ctrl and Exp mice
datadir_pre_exp = [xlsdir '\Pre_Surgery_Performance_Exp'];
datadir_pre_ctrl = [xlsdir '\Pre_Surgery_Performance_Ctrl'];

choosecb('MS_WM_EXP_cellbase');
datadir_post_exp = getpref('cellbase','datapath');

choosecb('MS_WM_CTRL_cellbase');
datadir_post_ctrl = getpref('cellbase','datapath');

% Extract behavioral data from Excel file
[mouse_pre_exp, mouse_post_exp, mouse_pre_ctrl, mouse_post_ctrl,...
   NWM7_pre, NWM9_pre, NWM15_pre, NWM7_post, NWM9_post, NWM15_post,...
   NWM19_pre, NWM21_pre, NWM22_pre, NWM19_post, NWM21_post, NWM22_post,...
   NWM19_fix, NWM21_fix, NWM22_fix] = extract_excel_data(xlsdir);

% Supplementary Fig S2: Pre vs Post surgery performance (Exp)
mouse_ids = {'NWM5','NWM7','NWM9','NWM15'};
pre_vs_post_perf(datadir_pre_exp, datadir_post_exp, 'exp', mouse_ids, resdir);

% Supplementary Fig S4: Pre vs Post surgery performance (Ctrl)
mouse_ids = {'NWM19', 'NWM21', 'NWM22'};
pre_vs_post_perf(datadir_pre_ctrl, datadir_post_ctrl, 'ctrl', mouse_ids, resdir);

% Supplementary Fig S1: Behavioral trends (Exp)
mouse_ids = {'NWM5','NWM7','NWM9','NWM15'};

figure('Units', 'centimeters', 'Position', [0, 0, 17, 25]); 

% Colors for each mouse - correct and incorrect
mouse_colors_correct = linspace(1, 0.5, size(mouse_ids, 2))' * [0.4, 1, 0.4]; % Green shades
mouse_colors_incorrect = linspace(1, 0.5, size(mouse_ids, 2))' * [1, 0.4, 0.4]; % Red shades

% S1.A - Pre surgery Performance of Correct and Incorrect Trials
height = 0.08;
subplot('Position',[0.08,0.89,0.669999999999998, height]); % Position 1: first plot
legend1_pos = [0.761534473774283,0.88,0.237947118134002,0.1];
hold on;
plotMouseData(mouse_pre_exp, mouse_colors_correct, mouse_colors_incorrect, mouse_ids, legend1_pos);
set(gca, 'TickLength', [0.015, 0.015]);
hold off;

% S1.B - Post surgery Performance of Correct and Incorrect Trials
subplot('Position',[0.08,0.72,0.67, height]); 
hold on;
legend2_pos = [0.761534473774283,0.71,0.237947118134002,0.1];
plotMouseData(mouse_post_exp, mouse_colors_correct, mouse_colors_incorrect, mouse_ids, legend2_pos);
set(gca, 'TickLength', [0.015, 0.015]);
hold off;

% S1.C - Stacked bar for NWM7, session perfoemance
subplot('Position',[0.08,0.55,0.669999999999998,height]); % Position 3: third plot
legend3_pos = [0.761534473774283,0.56,0.237947118134002,0.07];
plotBar(NWM7_pre, NWM7_post, legend3_pos,'Exp');
set(gca, 'TickLength', [0.015, 0.015]);

% Stacked bar for NWM9, session performance
subplot('Position',[0.08,0.38,0.669999999999998,height]);
legend4_pos = [0.761534473774283,0.39,0.237947118134002,0.07];
plotBar(NWM9_pre, NWM9_post, legend4_pos,'Exp');
set(gca, 'TickLength', [0.015, 0.015]);

% Stacked bar for NWM15, session performance
subplot('Position',[0.08,0.21,0.669999999999998,height]); % Position 5: fourth plot
legend5_pos = [0.761534473774283,0.22,0.237947118134002,0.07];
plotBar(NWM15_pre, NWM15_post, legend5_pos,'Exp');
set(gca, 'TickLength', [0.015, 0.015]);

% Empty subplot for Illustrator
subplot('Position',[0.08,0.04,0.28,height]); 
axis off; 

set(gca, 'TickLength', [0.015, 0.015]);
% Set font size for all text
set(findall(gcf, '-property', 'fontsize'), 'FontSize', 8);

% Save figure
saveas(gcf, fullfile(resdir, 'S1exp.svg'));
saveas(gcf, fullfile(resdir, 'S1exp.jpg'));

% Supplementary Fig S3: Behavioral trends (Ctrl)

figure('Units', 'centimeters', 'Position', [0, 0, 17, 25]); % Adjust figure size for vertical layout
mouse_ids = {'NMW19', 'NWM21', 'NWM22'};

% S3.A Pre surgery Performance of Correct and Incorrect Trials

% Colors for each mouse - correct and incorrect
mouse_colors_correct = linspace(1, 0.5, size(mouse_ids, 2))' * [0.4, 1, 0.4]; % Green shades
mouse_colors_incorrect = linspace(1, 0.5, size(mouse_ids, 2))' * [1, 0.4, 0.4]; % Red shades

height = 0.08;
subplot('Position',[0.08,0.89,0.669999999999998, height]); % Position 1: first plot
hold on;
legend1_pos = [0.761534473774283,0.88,0.237947118134002,0.1];
plotMouseData(mouse_pre_ctrl, mouse_colors_correct, mouse_colors_incorrect, mouse_ids, legend1_pos);
set(gca, 'TickLength', [0.015, 0.015]);
hold off;

% S3.B Post surgery Performance of Correct and Incorrect Trials

subplot('Position',[0.08,0.72,0.67, height]); 
hold on;
legend2_pos = [0.761534473774283,0.71,0.237947118134002,0.1];
plotMouseData(mouse_post_ctrl, mouse_colors_correct, mouse_colors_incorrect, mouse_ids, legend2_pos);
set(gca, 'TickLength', [0.015, 0.015]);
hold off;

% S3.C - Stacked bar for NWM19, session performance
subplot('Position',[0.08,0.55,0.669999999999998,height]); % Position 3: third plot
legend3_pos = [0.761534473774283,0.56,0.237947118134002,0.07];
plotBar(NWM19_pre, NWM19_post, legend3_pos,'Ctrl');
set(gca, 'TickLength', [0.015, 0.015]);

% Stacked bar for NWM21, session performance
subplot('Position',[0.08,0.38,0.669999999999998,height]);
legend4_pos = [0.761534473774283,0.39,0.237947118134002,0.07];
plotBar(NWM21_pre, NWM21_post, legend4_pos,'Ctrl');
set(gca, 'TickLength', [0.015, 0.015]);

% Stacked bar for NWM15, session performance
subplot('Position',[0.08,0.21,0.669999999999998,height]); % Position 4: fourth plot
legend5_pos = [0.761534473774283,0.22,0.237947118134002,0.07];
plotBar(NWM22_pre, NWM22_post, legend5_pos,'Ctrl');
set(gca, 'TickLength', [0.015, 0.015]);

% S3.D - Increasing fixation times: bar graph
subplot('Position',[0.08,0.04,0.28,height]); 
plot_fixationtime(NWM19_fix, NWM21_fix, NWM22_fix);
hold on;
set(gca, 'TickLength', [0.022, 0.015]);
hold off;

% Increasing fixation times: learning curve
subplot('Position',[0.46,0.04,0.28,height]);
legend6_pos = [0.7,0.07,0.237947118134002,0.06];
plot_fixationbar(NWM19_fix, NWM21_fix, NWM22_fix, legend6_pos);
hold on;
set(gca, 'TickLength', [0.022, 0.022]);
hold off;
set( findall(gcf, '-property', 'fontsize'), 'FontSize', 8);

% Save figure
saveas(gcf,fullfile(resdir, 'S3ctrl.svg'));
saveas(gcf,fullfile(resdir, 'S3ctrl.jpg'));

end 

function plotMouseData(mouse_data, colors_correct, colors_incorrect, mouse_ids, legend_pos)
% Plot performance curves (completed trials vs sessions)

    legendEntries = {};
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

        % Plot the curves based on the specified plot type
        plot(training_days, correct_trials, 'Color', colors_correct(mouse, :), 'LineWidth', 2);
                legendEntries{end + 1} = sprintf('%s Correct trials', mouse_ids{mouse});

        plot(training_days, incorrect_trials, 'Color', colors_incorrect(mouse, :), 'LineWidth', 2);
                legendEntries{end + 1} = sprintf('%s Error trials', mouse_ids{mouse});

    end
    
    % Labels and formatting
    xlabel('Training Session');
    ylabel('Completed Trials (%)','FontSize',10);
    ylim([0, 1]);
    yticks([0, 1]);
    legend(legendEntries, 'Position',legend_pos, 'NumColumns', 1,'Box', 'off','FontSize',8);
    box off;
    set(gca, 'TickDir', 'out');
end

function plotBar(mousePerfData_pre, mousePerfData_post, legend_pos, type)
% Plot stack bar session performance graphs

    if strcmp(type, 'Ctrl')
        mousePerfData = [mousePerfData_pre(1,:) mousePerfData_post(1,:);  mousePerfData_pre(2,:) mousePerfData_post(2,:);  mousePerfData_pre(3,:) mousePerfData_post(3,:);  mousePerfData_pre(4,:) mousePerfData_post(4,:)];
    else 
        correct_trials = sum(mousePerfData_pre(1:2, :), 1);
        incorrect_trials = sum(mousePerfData_pre(3:4, :), 1);
        aborted_trials = sum(mousePerfData_pre(5:6, :), 1);
        missed_trials = sum(mousePerfData_pre(7:8, :), 1);
    
        mousePerfData = [correct_trials mousePerfData_post(1,:);  incorrect_trials mousePerfData_post(2,:);  aborted_trials mousePerfData_post(3,:);  missed_trials mousePerfData_post(4,:)];
    end
    mousePerfData = vertcat(mousePerfData(1,:), mousePerfData(2,:), mousePerfData(3,:), mousePerfData(4,:));
    ba1 = bar(mousePerfData', 0.6, 'stacked');
    ylim([0 1]);
    ba1(1).FaceColor = [0 0.8 0]; % Correct
    ba1(2).FaceColor = [0.8 0 0]; % Incorrect
    ba1(3).FaceColor = [0.6 0.6 0.6]; % Abort
    ba1(4).FaceColor = [0.9 0.9 0]; % Timed Out
    ax=gca;
    surgery_time=size(mousePerfData_pre,2)+0.5;
    axis_position = get(ax, 'Position'); % [x, y, width, height] in normalized figure units
    
    x_norm_in_axis = (surgery_time - ax.XLim(1)) / (ax.XLim(2) - ax.XLim(1)); % X position in axis
    y_norm_in_axis = 1.15; % Position the arrow at the top of the axis (1 in normalized y-coordinates)
    x_fig = axis_position(1) + x_norm_in_axis * axis_position(3); % Convert x to figure space
    y_fig_start = axis_position(2) + y_norm_in_axis * axis_position(4); % Start y at top of axis
    y_fig_end = y_fig_start - 0.016; % Small downward offset for the arrow 
    annotation('textarrow', [x_fig, x_fig], [y_fig_start, y_fig_end], 'HeadStyle', 'vback1', 'LineWidth', 2, ...
        'Color', 'red', 'String','Surgery','FontSize',10);
    xlabel('Session');
    ylabel('Session performance (%)');
    yticks([0,1]);
    box off;
    set(gca, 'TickDir', 'out');
    legend('Correct trials', 'Error trials', 'Aborted trials', 'Missed trials', 'Position',legend_pos,'Box','off');
     
end

function pre_vs_post_perf(datadir_pre, datadir_post, txt, mice_names, resdir)
% Main wrapper for pre/post surgery comparison figures
    figure_width = 19;
    figure_height = 17; % in cm
    figure('Units', 'centimeters', 'Position', [2.09020833333333,6.40291666666667,figure_width,figure_height]);
    
    % Constants for layout
    vertical_gap = 0.1; % Vertical space between subplots
    margin_bottom = 0.08; % Bottom margin (to leave space for Illustrator edits)
    legend_spacing = 0.04;
    legend1_height = 0.12;
    legend2_height = 0.06;
    legend3_height = 0.06;
    legend1_width = 0.33; % Widest legend (4 lines)
    legend2_width = 0.3; % Narrower legend (2 lines)
    legend3_width = 0.3;
    % Calculate positions for subplots
    legend1_y = 1 - legend1_height - legend_spacing; % Top-aligned
    num_rows = 2;  % Number of subplot rows
    figure_height = legend1_y  - margin_bottom - vertical_gap - 0.04; % Total height for rows 2 & 3
    row_height = figure_height / num_rows; % Each subplot row gets half the height
    subplot_width = 0.25;
    row3_y = margin_bottom;
    row2_y = legend1_y - row_height - 0.015;
    % Legend positions (Row 1)
    legend1_x = 0.025;
    legend2_x = legend1_x + legend1_width + 0.025;
    legend2_y = 1 - legend2_height - legend_spacing; % Aligned to bottom of legend1
    legend3_x = 0.7;
    legend3_y = legend2_y; % Same height as legend2
    
    % Subplot positions (Rows 2 & 3)
    subplot1_x = 0.08; % Align to widest legend
    subplot2_x = subplot1_x + subplot_width + 0.08; % Align with legend 2
    subplot3_x = subplot2_x + subplot_width + 0.08; % Align with legend 3

    % Fig S2 or S4 A and B
    % A - Temporal delta between end of the mandatory central port fixation cue and the mouse
    % disengagement from it, presented as probability density and cumulative distribution
    % function.
    % B - Pre vs post surgery comparison of withdrawl percentage and
    % successfull/unsuccessful distribution on total trial.
    % For all figures: correct, uncorrect and withdrawed trials are
    % respectively green, red and gray, with light and dark hue for the before
    % and after surgey.
    [GvaluesBefore, BvaluesBefore, WvaluesBefore] = performance_extrapolation(datadir_pre,mice_names,'pre');
    [GvaluesAfter, BvaluesAfter, WvaluesAfter] = performance_extrapolation(datadir_post,mice_names,'post');
    
    % (In)Completed Trials ksdensity
    subplot('Position', [subplot1_x, row2_y, subplot_width, row_height]);
    legend1_pos = [legend1_x, legend1_y, legend1_width, legend1_height];
     before_vs_after_ksdensity({-WvaluesBefore,-WvaluesAfter,[GvaluesBefore  BvaluesBefore], [GvaluesAfter  BvaluesAfter] },...
        {[0.7, 0.7, 0.7],[0.5, 0.5, 0.5],[0.2, 0.8, 0.8],[0, 0.5, 1]},'a', legend1_pos, txt);
    
     % Correct trials ksdensity
    subplot('Position', [subplot2_x, row2_y, subplot_width, row_height]);
    legend2_pos = [legend2_x, legend2_y, legend2_width, legend2_height];
    before_vs_after_ksdensity({GvaluesBefore, GvaluesAfter}, {[0, 1, 0],[0, 0.5, 0]},'b', legend2_pos, txt)
    
    % Incorrect trials ksdensity
    subplot('Position', [subplot3_x, row2_y, subplot_width, row_height]);
    legend3_pos = [legend3_x, legend3_y, legend3_width, legend3_height];
    before_vs_after_ksdensity({BvaluesBefore, BvaluesAfter},{[1, 0, 0],[0.5, 0, 0]},'c', legend3_pos, txt)
    
    % Aborted
    subplot('Position', [subplot1_x, row3_y, subplot_width, row_height]);
    colours = {[0.7 0.7 0.7],[0.5 0.5 0.5]};
    before_vs_after_bar(colours,'a', WvaluesBefore, GvaluesBefore, WvaluesAfter, GvaluesAfter, BvaluesBefore, BvaluesAfter,txt)
    
    % Successfull trials
    subplot('Position', [subplot2_x, row3_y, subplot_width, row_height]);
    colours = {[0 1 0],[0 0.7 0]};
    before_vs_after_bar(colours, 'b', GvaluesBefore, BvaluesBefore, GvaluesAfter, BvaluesAfter, txt)
    
    % Unsuccessfull trials
    subplot('Position', [subplot3_x, row3_y, subplot_width, row_height]);
    colours = {[1 0 0],[0.7 0 0]};
    before_vs_after_bar(colours, 'c', BvaluesBefore, GvaluesBefore, BvaluesAfter, GvaluesAfter, txt)
    set( findall(gcf, '-property', 'fontsize'), 'fontsize', 10,'LineWidth', 1);
    
    if ~isempty(resdir)
        if strcmp(txt, 'ctrl')
            name = 'S4';
        else 
            name = 'S2';
        end
        fnm = fullfile(resdir, [name txt '.png']); % saving
        fnmm = fullfile(resdir, [name txt '.svg']); % saving
        saveas(gcf,fnm);
        saveas(gcf,fnmm);
    end
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

function before_vs_after_ksdensity(values_struct, line_colour,str, lgd_pos, txt)
% Plots kernel density estimates : ksdensity distribution plot of 2 events time distributions
% (values_struct) before and after surgery, represented with distinct colours (line_colour)

hold on
for i = 1 : length(values_struct)
    [f1,x1] = ksdensity(values_struct{i}); %
    plot(x1,f1,'color',line_colour{i},'LineWidth',2)
end
xlabel(sprintf('Time from central port disengage (s)'),'FontSize',8);
ylabel('Probability','FontSize',8.5);
if strcmp(str,'a')
    legend('Aborted trials before surgery','Aborted trials after surgery','Completed trials before surgery','Completed trials after surgery','Position',lgd_pos,'Box','off','FontSize',8);
    if strcmp(txt, 'exp')
        xlim([-2.5, 4]);
        xticks([-2.5,-2,0,2,4]);
        xticklabels({'','-2','0','2','4'});
        yticks([0,2,4]);
        ylim([0,4]);
    else 
        xlim([-1.5, 3]);
        xticks([-1.5,-1,0,1,2,3]);
        xticklabels({'','-1','0','1','2','3'});
        yticks([0,5,10,15]);
        ylim([0,15]);
    end
elseif strcmp(str,'b')
    legend('Correct trials before surgery','Correct trials after surgery','Position',lgd_pos,'Box','off','FontSize',8);
    if strcmp(txt, 'exp')
        xlim([-1, 4]);
        xticks([-1,0,2,4]);
        xticklabels({'', '0','2','4'});
        yticks([0,1,2,2.5]);
        yticklabels({'0', '1','2',''});
        ylim([0,2.5]);
    else 
        xlim([-1, 3]);
        xticks([-1,0,1,2,3]);
        yticks([0,3,5]);
        yticklabels({'0', '3','5'});
        ylim([0,5]);
    end
else 
    legend('Error trials before surgery','Error trials after surgery','Position',lgd_pos,'Box','off','FontSize',8);
    if strcmp(txt, 'exp')   
        xlim([-1, 4]);
        xticks([-1,0,2,4]);
        xticklabels({'', '0','2','4'});
        yticks([0,1,2,2.5]);
        yticklabels({'0', '1','2',''});
        ylim([0,2.5]);
    else 
        xticks([-1,0,1,2,3,4]);
        xlim([-1,4]);
        yticks([0,3,5]);
        ylim([0,5]);
        yticklabels({'0', '3','5'});
    end
end 
set(gca,'TickDir','out');
hold off
end

function before_vs_after_bar(bar_colour, str, MainBefore, SideBefore, MainAfter, SideAfter, SideBefore2, SideAfter2, txt)
% Plots bar graphs with error bars to compare specific performance outcomes before and after surgery  

    switch nargin % switching between percentage referred to completed trials only or all trials
        case 7 % Completed trials
            total_before = numel(MainBefore) + numel(SideBefore);
            total_after = numel(MainAfter) + numel(SideAfter);
            bef = numel(MainBefore) / total_before;
            aft = numel(MainAfter) / total_after;
        case 9 % Alls Trials
            total_before = numel(MainBefore) + numel(SideBefore) + numel(SideBefore2);
            total_after = numel(MainAfter) + numel(SideAfter) + numel(SideAfter2);
            bef = numel(MainBefore) / total_before;
            aft = numel(MainAfter) / total_after;
    end
    
    b_bar = [bef*100 aft*100];
    hold on
    bar(0,b_bar(1), 1, 'FaceColor',[1 1 1 ],'EdgeColor', bar_colour{1},'LineWidth',2);
    bar(3,b_bar(2), 1,'FaceColor', [ 1 1 1],'EdgeColor', bar_colour{2},'LineWidth',2); % second column color fix
    xticks([0,3]);
    if strcmp(str,'a')
        ylabel('All Trials (%)','FontSize',9,'VerticalAlignment','bottom', 'Units','centimeters','Position',[-0.56391527146962,2.632606677313645,0]);  
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
        ylabel('Completed Trials (%)','FontSize',9,'VerticalAlignment','bottom', 'Units','centimeters','Position',[-0.56391527146962,2.632606677313645,0]);
        ylim([0,90]);
        xticklabels({''});
        text(0,-12,{'Correct trials'; 'before surgery'},'HorizontalAlignment','center','FontSize',8);
        text(3,-12,{'Correct trials'; 'after surgery'},'HorizontalAlignment','center','FontSize',8);
        yticks([0,45,90]);
        yticklabels({'0','','90'});
    else 
        ylabel('Completed Trials (%)','FontSize',9,'VerticalAlignment','bottom', 'Units','centimeters','Position',[-0.56391527146962,2.632606677313645,0]);
        ylim([0,90]);
        xticklabels({''});
        text(0,-12,{'Error trials'; 'before surgery'},'HorizontalAlignment','center','FontSize',8);
        text(3,-12,{'Error trials'; 'after surgery'},'HorizontalAlignment','center','FontSize',8);
        yticks([0,45,90]);
        yticklabels({'0','','90'});
    end
    set(gca,'TickDir','out');
    hold off
end

function plot_fixationtime(NWM19, NWM21, NWM22)
% Line plot of fixation time trends

    % Fixation increase for Ctrl
    % Colors for mouse
    color19 = [0.678, 0.847, 0.902];
    color21 = [0, 0, 1];
    color22 = [0, 0, 0.545];
    hold on
    plot(NWM19,'color',color19,'LineWidth',2)
    plot(NWM21,'color',color21,'LineWidth',2)
    plot(NWM22,'color',color22,'LineWidth',2)
    plot([0,length(NWM22)], [0,0],'k')
    xlim([0,length(NWM22)]);
    ylim([0,1]);
    yticks([-2,0,1,2]);
    ylabel('Fixation length (s)','FontSize',9);
    xlabel('Training Session');
    set(gca,'TickDir','out','TitleHorizontalAlignment','left','TickLength',[0.02,0.02]);
    hold off;
end

function plot_fixationbar(NWM19, NWM21, NWM22, lgd)
% Bar plot comparing first and last day fixation times
    
    % Colors for mouse
    color19 = [0.678, 0.847, 0.902];
    color21 = [0, 0, 1];
    color22 = [0, 0, 0.545];
    hold on
    barmeans = [mean([NWM19(1) NWM21(1) NWM22(1)]) , mean([NWM19(end) NWM21(end) NWM22(end) ])];
    bar([0,3],barmeans, 0.6, 'FaceColor','none','EdgeColor',[0 0 0.2],'LineWidth',2)
    plot([0,3],[NWM19(1) NWM19(end)],'color',color19,'LineWidth',2)
    plot([0,3],[NWM21(1) NWM21(end)],'color',color21,'LineWidth',2)
    plot([0,3],[NWM22(1) NWM22(end)],'color',color22,'LineWidth',2)
    plot([0,3], [0,0],'k')
    ylabel('Fixation length (s)','FontSize',9);
    xticks([0,3]);
    xticklabels({''});
    xlim([-1.2,4.2]);
    ylim([0,1.5]);
    yticks([0,1,2]);
    text(0, -0.4,{'First'; 'training session'},'HorizontalAlignment','center','FontSize',9);
    text(3,-0.4,{'Last'; 'training session'},'HorizontalAlignment','center','FontSize',9);
    ax=gca;
    ax.TickDir='out';
    legend({'','NWM19','NWM21','NWM22'},'Box','off','NumColumns',1, 'Position', lgd);
    set(gca,'TickDir','out','TickLength',[0.02,0.02]);
    hold off
end

function [mouse_pre_exp, mouse_post_exp, mouse_pre_ctrl, mouse_post_ctrl,...
   NWM7_pre, NWM9_pre, NWM15_pre, NWM7_post, NWM9_post, NWM15_post,...
   NWM19_pre, NWM21_pre, NWM22_pre, NWM19_post, NWM21_post, NWM22_post,...
   NWM19_fix, NWM21_fix, NWM22_fix] = extract_excel_data(xlsdir)
% Reads data from Excel sheets
    
    % Pre & post surgery perf data for exp animals
    mouse_pre_exp=cell(1, size(4, 2) );
    mouse_pre_exp{1} = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'B85:AE88');
    mouse_pre_exp{2} = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'B96:AE99');
    mouse_pre_exp{3} = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'B109:BX112');
    mouse_pre_exp{4} = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'B122:AJ125');
    
    NWM7_pre = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'B96:AE103');
    NWM9_pre = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'B109:BX116');
    NWM15_pre = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'B122:AJ129');
    
    mouse_post_exp=cell(1, size(4, 2) );
    mouse_post_exp{1} = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'CP19:EB22');
    mouse_post_exp{2} = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'CP26:EL29');
    mouse_post_exp{3} = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'CP34:EF37');
    mouse_post_exp{4} = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'CP43:FX46');
    
    NWM7_post = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'CP26:EL29');
    NWM9_post = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'CP34:EF37');
    NWM15_post = xlsread([xlsdir '\MiceTime.xlsx'], 1, 'CP43:FX46');
    
    % Pre & post surgery perf data for ctrl animals
    mouse_pre_ctrl=cell(1, size(3, 2) );
    mouse_pre_ctrl{1} = xlsread([xlsdir '\MiceTime.xlsx'], 2, 'B8:P9');
    mouse_pre_ctrl{2} = xlsread([xlsdir '\MiceTime.xlsx'], 2, 'B14:O15');
    mouse_pre_ctrl{3} = xlsread([xlsdir '\MiceTime.xlsx'], 2, 'B21:Q22');
    
    NWM19_pre = xlsread([xlsdir '\MiceTime.xlsx'], 2, 'B8:P11');
    NWM21_pre = xlsread([xlsdir '\MiceTime.xlsx'], 2, 'B14:O17');
    NWM22_pre = xlsread([xlsdir '\MiceTime.xlsx'], 2, 'B21:Q24');
    
    mouse_post_ctrl=cell(1, size(3, 2) );
    mouse_post_ctrl{1} = xlsread([xlsdir '\MiceTime.xlsx'], 2, 'B31:AZ32');
    mouse_post_ctrl{2} = xlsread([xlsdir '\MiceTime.xlsx'], 2, 'B37:AZ38');
    mouse_post_ctrl{3} = xlsread([xlsdir '\MiceTime.xlsx'], 2, 'B43:BA44');
    
    NWM19_post = xlsread([xlsdir '\MiceTime.xlsx'], 2, 'B31:AZ34');
    NWM21_post = xlsread([xlsdir '\MiceTime.xlsx'], 2, 'B37:AZ40');
    NWM22_post = xlsread([xlsdir '\MiceTime.xlsx'], 2, 'B43:BA46');
    
    % Fixation times - Ctrl mice 
    NWM19_fix = xlsread([xlsdir '\MiceTime.xlsx'], 4, 'B1:AB1');
    NWM21_fix = xlsread([xlsdir '\MiceTime.xlsx'], 4, 'B5:AD5');
    NWM22_fix = xlsread([xlsdir '\MiceTime.xlsx'], 4, 'B8:AE7');
end