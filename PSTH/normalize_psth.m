function [normalizedDataWithoutNaN, normalizedPSTH, numPSTH, delayDataWithoutNaN] = normalize_psth(spsth_data, time, baselinePeriod, delayPeriod)
%NORMALIZE_PSTH Normalize spike PSTH data using a baseline period.
%   [NORMALIZEDDATAWITHOUTNAN, NORMALIZEDPSTH, NUMPSTH, DELAYDATAWITHOUTNAN] = 
%   NORMALIZE_PSTH(SPSTH_DATA, TIME) normalizes the input spike PSTH data matrix
%   SPSTH_DATA using a baseline period of [-1.6, 0] seconds by default, and
%   returns the normalized data.
%
%   [___] = normalize_psth(SPSTH_DATA, TIME, BASELINEPERIOD, DELAYPERIOD) allows
%   specification of custom baseline and delay periods.
%
%   Inputs:
%       SPSTH_DATA       - Spike PSTH data matrix. Size: NUMPSTH x NUMTIMEPOINTS
%       TIME             - Time vector corresponding to the PSTH data. Size: 1 x NUMTIMEPOINTS
%       BASELINEPERIOD   - [start, end] (optional), baseline period in seconds.
%                          Default: [-1.6, 0]
%       DELAYPERIOD      - [start, end] (optional), delay period in seconds to extract.
%                          Default: [0.2, 1]
%
%   Outputs:
%       NORMALIZEDDATAWITHOUTNAN - Normalized PSTH data with NaN rows removed.
%                                  Size: NUMVALIDPSTH x NUMTIMEPOINTS
%       NORMALIZEDPSTH           - Full normalized PSTH data (including NaNs).
%                                  Size: NUMPSTH x NUMTIMEPOINTS
%       NUMPSTH                  - Number of PSTHs in the input data.
%       DELAYDATAWITHOUTNAN      - Extracted and cleaned delay period data.
%                                  Size: NUMVALIDPSTH x NUMDELAYPOINTS
%
%   Example:
%       normalizedData = normalize_psth(spsth_data, time);
%
%       % With custom baseline and delay periods:
%       normalizedData = normalize_psth(spsth_data, time, [-1, 0], [0.1, 0.8]);
%
%
%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025

    % Set default values for optional parameters
    if nargin < 3
        baselinePeriod = [-1.6, 0]; % Default baseline period
    end
    if nargin < 4
        delayPeriod = [0.2, 1]; % Default delay period
    end

    % Define baseline period for normalization
    baselineStart = baselinePeriod(1); % Start time of baseline period (seconds)
    baselineEnd = baselinePeriod(2); % End time of baseline period (seconds)
    baselineIndices = find(time >= baselineStart & time <= baselineEnd);

    % Check if baseline period indices are valid
    if isempty(baselineIndices)
        error('Baseline period not found in the provided time vector.');
    end

    % Get dimensions of spsth_data
    [numPSTH, numTimePoints] = size(spsth_data);
    normalizedPSTH = nan(numPSTH, numTimePoints);

    % Loop through each PSTH to normalize using baseline period
    for iPSTH = 1:numPSTH
        meanBaseline = mean(spsth_data(iPSTH, baselineIndices)); % Mean of baseline period
        stdBaseline = std(spsth_data(iPSTH, baselineIndices)); % Standard deviation of baseline period
        normalizedPSTH(iPSTH, :) = (spsth_data(iPSTH, :) - meanBaseline) / stdBaseline; % Z-score normalization
    end

    % Define delay period
    delayStart = delayPeriod(1); % Start time of delay period (seconds)
    delayEnd = delayPeriod(2); % End time of delay period (seconds)

    % Extract and clean delay data
    delayIndices = find(time >= delayStart & time <= delayEnd);

    % Check if delay period indices are valid
    if isempty(delayIndices)
        error('Delay period not found in the provided time vector.');
    end

    delayData = normalizedPSTH(:, delayIndices);
    delayDataWithoutNaN = delayData(sum(isnan(delayData), 2) == 0, :); % Remove rows with NaN values

    % Remove NaN rows from the normalized PSTH
    normalizedDataWithoutNaN = normalizedPSTH(sum(isnan(normalizedPSTH), 2) == 0, :);

end
