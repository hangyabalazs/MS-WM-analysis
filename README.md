# MATLAB code for analyzing medial septal single-unit activity during a working memory task

## Description

This repository contains MATLAB code for analyzing medial septal (MS) single-unit recordings collected from mice performing a delayed response two-alternative forced choice (DR-2AFC) working memory (WM) task and a control experiment. The full analysis pipeline can be executed by running the main function: `ms_wm_main.m`.

## Content

- **ACG/** – Functions for autocorrelogram (ACG) analysis
- **PSTH/** – Functions for peri-stimulus time histogram (PSTH) analysis
- **Behaviour/** – Functions and data for behavioral analysis
  - Contains behavioral datasets used in the study
- **Main script:** `ms_wm_main.m` – Executes full analysis pipeline
- **Required CellBases:**
  - `MS_WM_EXP_cellbase` – Neural data from mice performing the DR-2AFC WM task
  - `MS_WM_CTRL_cellbase` – Neural data from mice performing the DR-2AFC control task
- **Generated output:**
  - Results are saved in the folder `MS_WM_results/`
  - Intermediate files like `MS_WM_data_filtered.mat`, and `ROC_pvalues.mat` are created if not present

## Installation

1. Download all `.m` files from this repository.
2. Add the following directories (and their subdirectories) to your MATLAB path:
   - `ACG/`
   - `PSTH/`
   - `Behaviour/`
3. Download and install [CellBase](https://github.com/hangyabalazs/CellBase ).
4. Mount the required CellBases:
     - Run `mountcb.m` and name the mounted CellBases:
         - `MS_WM_EXP_cellbase` for the working memory experiment
         - `MS_WM_CTRL_cellbase` for the control experiment
     - When prompted to "locate CellBase database file", select the `Cleansed_CB_EXP.mat` or `Cleansed_CB_CTRL.mat` file from each dataset
     - Choose `'CellBase'` for the timestamp conversion option

> A typical installation takes about 5–10 minutes, depending on dataset size and system performance.

## Demo / Example Usage

To run the full analysis:

1. Ensure that the CellBases are properly mounted as described above.
2. Make sure the `Behaviour/` directory contains the necessary behavioral data.
3. In the MATLAB Command Window, run:
   ms_wm_main
The results will be saved in a folder named MS_WM_results/ inside your current working directory directory.
If intermediate files like MS_WM_data_filtered.mat, and ROC_pvalues.mat already exist in the appropriate directories, the analysis will run faster.
 

## Instructions for Use

- Replace the behavioral data in the Behaviour/ folder with your own
- Update the mounted CellBase names in ms_wm_main.m to match your experimental datasets
- All figures and statistical outputs are automatically saved in the MS_WM_results/ directory

## Dependencies
- CellBase : Required for handling neural data - https://github.com/hangyabalazs/CellBase
- hangya-matlab-code package - https://github.com/hangyabalazs/Hangya-Matlab-code
- MATLAB Statistics and Machine Learning Toolbox


## System Requirements
- The code was developed and tested on the following system:

     - Operating System: Windows 10 Pro 64-bit
     - Processor: Intel(R) Core(TM) i7-6500U CPU @ 2.50GHz
     - RAM: 16 GB
     - MATLAB Version: R2023a
     - Toolboxes: Statistics and Machine Learning Toolbox


## Contact
For any questions, bug reports, or general feedback, please contact: hangya.balazs@koki.hun-ren.hu
