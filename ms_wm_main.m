function ms_wm_main
%MS_WM_MAIN Main analysis pipeline for MS WM experiment
%
% This function runs the complete analysis pipeline for the Medial Septum 
% Working Memory experiment. It processes behavioral data, performs PSTH 
% and ACG analyses, statistical testing, and generates all main figures 
% and supplementary figures.
%
% It sets up directory paths and creates results directory, generates 
% behavioral analysis plots (Fig. 1, Fig. S1-4), processes electrophysiological
% data or loads existing processed data, generates statistical analysis plots
% (Fig. 2), creates grouping analysis plots (Fig. 3, Fig. S5), performs ROC 
% analysis and generates Fig. 4 with Fig. S6-7, and analyzes theta rhythmicity
% using ACGs (Fig. 5, Fig. S8-10). It saves all main and supp. figures in
% main results directory, except for Fig. 3A-B, 3C & S5C which are saved in
% under \grouping folder.
%
% Example:
%   ms_wm_main()
%
% See also: Behaviour_panel, process_ms_wm_data, stat_panel_fig2, 
%           figure3_panel, roc_analysis, stat_panel_fig4, 
%           theta_rhythmicity_panel
%
%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025

    % Directories
    global datapath
    datapath = fileparts(which(mfilename));
    behav_data_fnm = fullfile(datapath,'Behaviour\'); % Behaviour data directory
    resdir = fullfile(datapath,'MS_WM2'); % Results directory
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
    
    % Figure 2 - Task phases stats
    stat_panel_fig2(resdir, cleaned_data);
     
    % Figure 3 - Delay response stat-based cell categories + S5 (ctrl)
    figure3_panel(resdir, [resdir '\grouping'], cleaned_data, 1); 
     
    % Fig 4 - All neurons & per category stats + S6-7
    rocdir = [resdir '\grouping\roc\submean']; % Directory for ROC results

    try
        load([rocdir '\ROC_pvalues.mat']); % Load results if available
    catch 
        warning('Cannot find processed data file. Running process again.');
        roc_analysis(resdir, [0.2 1], cleaned_data, 'submean'); % Run ROC analysis
        load([rocdir '\ROC_pvalues.mat']); % Load results
    end
    
    stat_panel_fig4(cleaned_data, ROC, ROCtime, pvalues, resdir, 'submean');
    
    % Figure 5 - Theta rhythmicity grouping + S8-10
    theta_rhythmicity_panel(resdir, cleaned_data, 'submean');

end