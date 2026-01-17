function ms_wm_main
%MS_WM_MAIN Main analysis pipeline for MS WM experiment
%
% MS_WM_MAIN runs the complete analysis pipeline for the Medial Septum 
% Working Memory experiment. It processes behavioral data, performs PSTH 
% and ACG analyses, statistical testing, and generates all main figures 
% and supplementary figures.
%
% It sets up directory paths and creates results directory, generates 
% behavioral analysis plots (Fig. 1, Fig. S1-2), processes electrophysiological
% data or loads existing processed data, generates statistical analysis plots
% (Fig. 3), creates grouping analysis plots (Fig. 4, Fig. S3), and clustering 
% analysis panels (Fig. S4-5), performs ROC  analysis and generates Fig. 5-6 
% with Fig. S6-7, and analyzes theta rhythmicity using ACGs (Fig. 7, Fig. S8-9). 
% It saves all main and supp. figures in main results directory, except for 
% Fig. 4A-B, 4C & S3C which are saved under \grouping folder and Fig.
% S5-6 which are saved under \clustering folder.
%
% Example:
%   ms_wm_main()
%
% See also: BEHAVIOUR_PANEL, PROCESS_MS_WM_DATA, STAT_PANEL_FIG3, 
%           FIGURE4_PANEL, ROC_ANALYSIS, STAT_PANEL_FIG5_6, 
%           THETA_RHYTHMICITY_PANEL

%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025

    % Directories
    global datapath
    datapath = fileparts(which(mfilename));
    behav_data_fnm = fullfile(datapath,'Behaviour\'); % Behaviour data directory
    resdir = fullfile(datapath,'MS_WM_results'); % Results directory
    if ~isfolder(resdir)
        mkdir(resdir); % Create directory if it doesn't exist
    end
    
    % Figure 1 - Behaviour plots
    Behaviour_panel(behav_data_fnm, resdir);
    
    % SPSTH analysis, ACG analysis, Delay response stat-based categorization
    try
        load([resdir '\MS_WM_data_filtered.mat']); % Load results if available
    catch 
        warning('Cannot find processed data file. Running process again.');
        cleaned_data = process_ms_wm_data(resdir); % Run analyses
    end
    
    % Figure 3 - Task phases stats
    stat_panel_fig3(resdir, cleaned_data);

    % Figure 4 - Delay response stat-based cell categories + S3 (ctrl)
    figure4_panel(resdir, [resdir '\grouping'], cleaned_data, 1); 

    % Clustering panels (S4-5)
    try
        load([resdir '\MS_WM_clustering_results.mat']); % Load clustering results if available
    catch 
        warning('Cannot find clustering data file. Locate file and run again.');
    end
    figure4_panel(resdir, [resdir '\clustering'], cleaned_data, 1, 'clusters_exp',...
    clusters_assign_exp, 'clusters_ctrl', clusters_assign_ctrl); 
 
    % Fig 5-6 - All neurons & per category stats + S6-7
    rocdir = [resdir '\grouping\roc\submean']; % Directory for ROC results
    try
        load([rocdir '\ROC_pvalues.mat']); % Load results if available
    catch 
        warning('Cannot find processed data file. Running process again.');
        roc_analysis(resdir, [0.2 1], cleaned_data, 'submean'); % Run ROC analysis
        load([rocdir '\ROC_pvalues.mat']); % Load results
    end
    stat_panel_fig5_6(cleaned_data, ROC, ROCtime, pvalues, resdir, 'submean');

    % Figure 7 - Theta rhythmicity grouping + S8-9
    theta_rhythmicity_panel(resdir, cleaned_data, 'submean');

end