process_ms_wm_data.m
create spsth data, acg (theta index) and categorize cells based on delay response
requires: ultimate_psth_wm.m
filterTrials_wm.m (includes filter for trials with short fixation time)
normalize_psth.m
acgmod.m
wm_responsesorter.m

Fig 1:
Behavior_panel.m

supp behavior figs: 
Behavior_panel_supp.m

needs folders containing bpod data for presurgery performance figs for Exp & ctrl, cellbases for postsurgery performance + excell file (MiceTime) (created by NS and modified
by MA, use latest version by MA updated with ctrl mice data).

Fig 2

Fig 2A-D produced outside of matlab
Fig 2E-J : stat_panel_fig2
requires processed data struct produced by process_ms_wm_data.m
barmeanstat.m
cdf_fig.m
spsth_parts

Fig 3
require speaker image in the root directory + processed psth data
grouping_panel_exp.m
grouping_panel_ctrl.m
clustering_panel_exp.m
clustering_panel_ctrl.m
rquire viewcell2b
ultimate_psth
chiSquareTest.m