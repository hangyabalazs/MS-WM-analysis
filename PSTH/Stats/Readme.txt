stat_panel_fig3 is called in the main function that creates statistics panel in Fig2E-J. It requires spsth_parts, stat_categories.m, barmeanstat.m, & cdf_fig.m 
  -spsth_parts calculates mean FR during delivery feedback to be used in Fig2I-J. This can take a while.
  -stat_categories.m plots mean FR comparison between Exp and Ctrl across task phases per response category. (Fig S4)
  -barmeanstat plots bar graphs for mean FR.
  -cdf_fig calculates and plots cumulative density functions.

stat_panel_fig5_6 is called in the main function that creates statistics panel in Fig5 & Fig6. It requires barmeanstat.m, cdf_fig.m, meanbargroups.m, meanSE_bar.m & roc_analysis.m
  -meanbargroups plots bar graphs of mean FR per group requires meanSE_bar.
  -roc_analysis runs ROC analysis on data. requires categorization matrices. results in Fig6E-F & S5.