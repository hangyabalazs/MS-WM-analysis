process_ms_wm_data.m performs psth analysis aligned to cue onset (FixationBeginning) using ultimate_psth_wm.m (requires filterTrials_wm and psth_stats_wm),
normalizes psth data (normalize_psth.m) according to baseline [-1.6 0], performs sorting of neurons based on statistically significant activity during the delay [0.2 1] 
and saves results in CellBase under property 'delay_response', and performs acg analysis to obtain theta index for each neuron. 
Saved results are cleaned_data, a struct containing normalzied spsth data, delay response categorization matrices, theta indices, psth time vector and cellids. 
This will be used in subsequent function calls.

figure3_panel is the main function used in ms_wm_main. 
It plots the main Figure 3A-B (spsth heatmap, averages per group & piechart) for stat based grouping.
It requires cleaned_data issued from process_ms_wm_data which has normalized spsth data & grouping/clustering matrices.
It calls grouping_panel_exp & grouping_panel_ctrl to plot Fig 3C and S5C for heatmap per group, average per group,
raster & psth plots for example cells from groups. The creation of these panels can cause matlab to crash which is minimized by having these functions separate.
It also performs chi square tests and saves the results to an excel file (StatResults), therefore requires chiSquareTest.m.
