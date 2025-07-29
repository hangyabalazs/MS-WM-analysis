Behaviour_panel.m - Main function calls excel data (

MiceTime) for plotting figure 1 (behaviour data for exp mice).
It includes:
- Task timeline diagram
- Temporal delta between fixation end and tone offset
- Performance over training sessions (pre- and post-surgery)

Behaviour_panel_supp.m - produces supplementary figures (e.g., S1–S4) for exp & ctrl mice that expand on the main behavioral analysis, uses same excel data file. 
It includes:
- Line plots of correct vs incorrect trials over training
- Stacked bar charts showing session performance
- Statistical comparisons before and after surgery
- Fixation duration trends for control animals

Inputs:
The functions expect:
- 'xlsdir': directory where the Excel file named `MiceTime.xlsx` is located.
- 'resdir': results directory.