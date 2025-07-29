function mean_FR = spsth_parts(datapath)
%SPSTH_PARTS computes mean FR for Reward & Punishment for Exp and Ctrl.
%  SPSTH_PARTS(DATAPATH) computes the mean firing rate (FR) for neurons in 
%  both Reward and Punishment alignment conditions across Experimental (Exp) 
%  and Control (Ctrl) datasets.
%
% Inputs: 
%     datapath - Full path to the directory where results will be saved or loaded from.
%
% Outputs: 
%     mean_FR - A 2x2 cell array containing mean firing rates for each
%     neuron under different conditions.
%         mean_FR{1,1}: Experimental dataset, Reward-aligned
%         mean_FR{1,2}: Experimental dataset, Punishment-aligned
%         mean_FR{2,1}: Control dataset, Reward-aligned
%         mean_FR{2,2}: Control dataset, Punishment-aligned
% Example: 
%      mean_FR = spsth_parts('/path/to/experiment/results/');
%
% Notes:
%     The function checks for the presence of a rasterpsth_.mat file in subdirectories named after the condition (e.g., RewardExp, PunishmentCtrl).
%     If results are missing, the user is prompted to run the analysis.
%
% See also : ULTIMATE_PSTH_WM, CHOOSECB, LOADCB, SELECTCELL.
%
%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025

    % Set up dataset names and alignment conditions
    datasets = {'Exp', 'Ctrl'};
    alignments = {'Reward', 'Punishment'};
    
    mean_FR = cell(2,2);

    % Loop through each dataset (Experimental and Control)
    for d = 1:length(datasets)
        if strcmp(datasets{d}, 'Exp')
            cellbase_name = 'MS_WM_EXP_cellbase';
        else
            cellbase_name = 'MS_WM_CTRL_cellbase';
        end
        choosecb(cellbase_name);
        loadcb;  % Reload database
        cellids = selectcell('"ID">20&"Lratio"<0.15&ismember("Anatomy",''MS'')');

        % Loop through Reward and Punishment alignments
        for a = 1:length(alignments)
            resdir1 = fullfile(datapath, [alignments{a}, datasets{d}]);
            
            % Check if directory and results file exist
            datafile = fullfile(resdir1, 'rasterpsth_.mat');
            if exist(datafile, 'file')
                fprintf('Results already exist for %s-%s. Skipping analysis.\n', datasets{d}, alignments{a});
            else
                if ~isfolder(resdir1)
                    mkdir(resdir1);
                end
                
                % Ask the user if they want to rerun the analysis
                userResponse = input(['No results found for ', datasets{d}, '-', alignments{a}, '. Run analysis? (y/n): '], 's');
                if strcmpi(userResponse, 'y')
                    % Run ultimate_psthloop
                    fprintf('Running analysis for %s-%s...\n', datasets{d}, alignments{a});
                    ultimate_psth_wm(cellids, 'trial', [alignments{a}, 'Beginning'], [-2 6], ...
                        'parts', '#Outcome', 'resdir', resdir1, 'dt', 0.001, 'sigma', 0.02, ...
                        'isadaptive', 0, 'maxtrialno', Inf, 'relative_threshold', 0.01, ...
                        'shevent', {'RewardBeginning', 'PunishmentBeginning', 'TrialEnd'},...
                        'display', false, 'event_filter', 'lowfixation_wm');
                else
                    fprintf('Skipping analysis for %s-%s.\n', datasets{d}, alignments{a});
                    continue;
                end
            end

            % Load and process results
            load(datafile, 'spsth_all', 'time', 'EventTimes');
            mean_FR{d,a} = compute_mean_FR(spsth_all, time, EventTimes, lower(alignments{a}));
            
            % Save computed mean firing rate
            fprintf('Mean firing rate saved for %s-%s.\n', datasets{d}, alignments{a});
        end
    end
    save(fullfile(datapath, 'mean_FR.mat'), 'mean_FR');
end

function mean_FR = compute_mean_FR(spsth, time, EventTimes, type)

    numCells = size(spsth, 1);
    mean_FR = nan(numCells, 1); % Preallocate for efficiency

    for iC = 1:numCells
        trialEnds = EventTimes(iC).TrialEnds;
        
         % Check if trialEnds is all NaN
        if all(isnan(trialEnds))
            endTime = 1; % Default to 1s if all values are NaN
        else
            if strcmp(type, 'reward')
                endTime = median(trialEnds, 'omitnan'); % Ignore NaNs in median calculation
            elseif strcmp(type, 'punishment')
                if median(trialEnds, 'omitnan') > 2
                    endTime = 1; % If TrialEnds > 2s, use 1s
                else
                    endTime = median(trialEnds, 'omitnan');
                end
            end
        end
        
        % Find the index in time that corresponds to endTime
        idx_start = find(time == 0, 1);
        idx_end = find(time >= endTime, 1);

        % Compute mean firing rate up to endTime
        if ~isempty(idx_end)
            if strcmp(type, 'reward')
                try
                    mean_FR(iC) = mean(spsth{iC, 1}(idx_start:idx_end)); 
                catch
                    continue;
                end
            else 
                try
                    mean_FR(iC) = mean(spsth{iC, 2}(idx_start:idx_end)); 
                catch
                    continue;
                end
            end
        end
    end
end
