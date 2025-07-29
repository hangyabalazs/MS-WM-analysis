function roc_analysis(resdir, testwin, cleaned_data, analysis_type)
%ROC_ANALYSIS Performs ROC analysis on MS WM data by cell response category.
%   ROC_ANALYSIS(RESDIR, TESTWIN, CLEANED_DATA, ANALYSIS_TYPE) performs
%   ROC analysis on experimental and control data grouped by delay response
%   categories. It calculates ROC values over time windows and saves the results
%   in a mat file.
%
%   Inputs:
%       RESDIR         - Directory to save output figures and data
%       TESTWIN        - 1x2 vector defining the time window for ROC
%       analysis (delay window).
%       CLEANED_DATA    - Struct containing preprocessed data:
%                         - delay_response.ResponseCategoriCtrl
%                         - delay_response.ResponseCategoriExp
%                         - cellids.cellids_exp_cl
%                         - cellids.cellids_ctrl_cl
%       ANALYSIS_TYPE   - Type of preprocessing applied to data:
%                         'submean' : subtract baseline mean
%                         'norm'    : subtract and divide by baseline std
%                         'none'    : raw values
%
%
% See also: ULTIMATE_PSTH.
%
%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025

    close all
    rocdir=[resdir '\grouping\roc\' analysis_type];
    if ~isfolder(rocdir)
        mkdir(rocdir)
    end

    % Extract variables 
    ResponseCategoriCtrl = cleaned_data.delay_response.ResponseCategoriCtrl;
    ResponseCategoriExp = cleaned_data.delay_response.ResponseCategoriExp;
    cellids_exp = cleaned_data.cellids.cellids_exp_cl;
    cellids_ctrl = cleaned_data.cellids.cellids_ctrl_cl;

    % Clean categories matrices
    MC = categorical(ResponseCategoriCtrl);
    ME = categorical(ResponseCategoriExp);

    % Define parameters
    time=linspace(-2,6,8001);
    baseline_start = -1.6;
    baseline_end = 0;
    dt = 0.001;
    wns = 0.05 / dt;   % 50 ms window (in data points)
    sht = 0.025 / dt;   % shift between overlapping windows (in data points)

    % Define number of windows
    nmw = floor((diff(testwin) / dt - wns) / sht) + 1;   % number of overlapping windows 

    % Separate cell IDs according to groups
    group_names = {'Inh', 'Act', 'Inh-Act', 'Act-Inh', 'NonResp'};
    Matrix = {MC, ME};

    idexp = cell(1, 5);
    idct = cell(1, 5);
    for g = 1:length(group_names)
        groupIndicesE = (ME == group_names{g});
        idexp{g} = cellids_exp(groupIndicesE, :);
        groupIndicesC = (MC == group_names{g});
        idct{g} = cellids_ctrl(groupIndicesC, :);
    end

    % ROC
    X = cell(1, 2);
    for m = 1:length(Matrix)
        if m == 1
            choosecb('MS_WM_CTRL_cellbase');
            cellids = idct;
        else
            choosecb('MS_WM_EXP_cellbase');
            cellids = idexp;
        end

        X{1, m} = cell(1, length(group_names));

        for g = 1:length(group_names)
            % Set the number of windows 
            X{1, m}{g} = nan(length(cellids{1, g}), nmw);

            for c = 1:length(cellids{1, g})
                try
                    [~, ~, ~, spt] = ultimate_psth_wm(cellids{1, g}{c, 1}, 'trial', 'FixationBeginning', [-2 6], 'dt', 0.001, 'sigma', 0.02, 'isadaptive',0, 'maxtrialno',Inf, 'relative_threshold',0.01,...
                       'event_filter', 'lowfixation_wm', 'display', false, 'issave', false);

                    % Find spt in delay window
                    spt_delay = spt(:, find(time >= testwin(1) & time <= testwin(2)));

                    baseline_bins = find(time >= baseline_start & time <= baseline_end);
                    spt_baseline = spt(:, baseline_bins);  % Extract baseline spikes

                    mean_baseline = mean(spt_baseline,2);  % Normalized baseline spike rate
                    std_baseline = std(spt_baseline, 0, 2);  % Normalized baseline spike rate

                    for w=1:nmw
                        inx1 = (w - 1) * sht + 1;
                        inx2 = inx1 + wns - 1;
                        % Compute data for the whole delay window
                        if strcmp(analysis_type,'submean')
                            X{1, m}{g}(c, w) = mean(mean(spt_delay(:, inx1:inx2),2) - mean_baseline);
                        elseif strcmp(analysis_type,'norm')
                            X{1, m}{g}(c, w) = mean(mean(spt_delay(:, inx1:inx2),2) - mean_baseline ./ std_baseline);
                        else 
                            X{1, m}{g}(c, w) = mean(mean(spt_delay(:, inx1:inx2),2));
                        end
                    end
                    disp(['Cell #' num2str(c) ' / ' num2str(length(cellids{1, g})) ' done......'])
                catch
                    disp('error');
                    continue;
                end
            end
        end
    end

    fnmm= ['ROC_data' analysis_type '.mat'];
    save(fullfile(rocdir,fnmm),'X','-mat');
    
    % ROC plot
    ROCtime = linspace(testwin(1), testwin(2),nmw);
    ROC=nan(5,nmw);
    SE=nan(5,nmw);
    pvalues=nan(5,nmw);
    for g = 1:length(group_names)
        for w = 1:nmw
            [ROC(g,w),pvalues(g,w),SE(g,w)]=rocarea(X{1,1}{1,g}(:,w),X{1,2}{1,g}(:,w),'bootstrap',5000); 
        end
        figure   % plot
        yyaxis left 
        b=bar(ROCtime,pvalues(g,:)); 
        b.FaceAlpha=0.5; 
        ylabel('p value');
        hold on; 
        yyaxis right 
        h=plot(ROCtime,ROC(g,:)); 
    
        ylabel('ROC','Rotation',-90);
        allChildren=get(gca, 'Children'); 
        newOrder=[setdiff(allChildren,h);h]; 
        set(gca, 'Children', newOrder)
        set(gca, 'TickDir', 'out', 'Box', 'off');
        title([group_names(g)]);
    
        fnm=[ 'ROC___' num2str(g) '.png'];
        saveas(gcf, fullfile(rocdir,fnm ));
        fnm2=[ 'ROC___' num2str(g) '.fig'];
        saveas(gcf, fullfile(rocdir,fnm2 ));
        close(gcf); 
    end
    save([rocdir '\ROC_pvalues.mat'], 'pvalues', 'ROC','ROCtime');

end
