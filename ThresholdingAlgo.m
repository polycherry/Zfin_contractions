function [signals, avgFilter, stdFilter] = ThresholdingAlgo(y, lag, threshold, influence)
% THRESHOLDINGALGO - Robust peak signal detection using dynamic thresholding
%
% This function implements a robust thresholding algorithm based on z-scores 
% and moving statistics to convert noisy signals into binary form (1 = movement, 
% 0 = no movement). The algorithm is particularly useful for detecting embryonic 
% contractions where signal amplitude varies significantly across developmental 
% stages.
%
% The algorithm constructs separate moving mean and standard deviation values 
% such that detected signals do not corrupt the threshold. This makes it robust 
% to varying signal amplitudes and noise levels.
%
% INPUTS:
%   y         - Input signal vector (1D array of signal values)
%               Typically squared frame-subtraction values from video analysis
%
%   lag       - Number of previous data points used to calculate moving statistics
%               Represents the window size for the moving mean and std deviation
%               Recommended: 30 (equivalent to 3 seconds at 10 fps)
%               Higher values = smoother detection but less responsive to changes
%
%   threshold - Z-score threshold for signal detection (number of standard deviations)
%               A point is considered a signal if it deviates by this many std devs
%               Recommended: 10 (based on validation with manual counting)
%               Higher values = fewer false positives but may miss weak signals
%
%   influence - Controls how detected signals affect future thresholds (0 to 1)
%               0 = detected signals don't influence threshold (most conservative)
%               1 = detected signals fully influence threshold (least conservative)
%               Recommended: 0.5 (balanced approach)
%               Also called "smoothing coefficient" in some contexts
%
% OUTPUTS:
%   signals   - Binary signal vector (same length as y)
%               1 = signal detected (movement/contraction occurring)
%               0 = no signal (no movement)
%
%   avgFilter - Moving average filter values over time
%               Useful for debugging and understanding threshold adaptation
%
%   stdFilter - Moving standard deviation filter values over time
%               Useful for debugging and understanding threshold adaptation
%
% ALGORITHM OVERVIEW:
%   1. Initialize filters using first (lag+1) data points
%   2. For each new data point:
%      a. Check if it deviates from moving average by threshold*std
%      b. If yes and above average: signal = 1
%      c. If yes but below average: signal = 0 (noise suppression)
%      d. If no: signal = 0
%   3. Update filtered signal with influence factor
%   4. Recalculate moving statistics using filtered signal
%
% VALIDATION:
%   This algorithm was validated against manual contraction counting using 
%   12 random embryos. The optimal parameters (lag=30, threshold=10, 
%   influence=0.5) resulted in minimal false positives/negatives with a 
%   sum of squared differences of 18.
%
% EXAMPLE:
%   % Detect peaks in a noisy signal
%   signal = [your_frame_subtraction_data];
%   [binary_signal, avg, std] = ThresholdingAlgo(signal, 30, 10, 0.5);
%   
%   % Find contraction events
%   [peaks, locs] = findpeaks(binary_signal, 'MinPeakDistance', 10);
%   num_contractions = length(peaks);
%
% REFERENCES:
%   Based on peak detection algorithms using moving statistics and z-scores.
%   Adapted for embryonic contraction detection with specific parameters
%   optimized through systematic validation.
%
% See also: FINDPEAKS, MEAN, STD
%
% Author: Vasishta Polisetty
% Developed under guidance of: Urvashi Jha (PhD student)
% Lab: Vatsala Thirumalai Lab, NCBS Bengaluru
% Date: 2017
% Updated: February 2026 for GitHub release

%% Initialize output arrays and filters
% Pre-allocate binary signal array for efficiency
signals = zeros(length(y), 1);

% Initialize the filtered signal with first (lag+1) points
% This forms the baseline for comparison
filteredY = y(1:lag+1);

% Calculate initial moving average and standard deviation
% These will be updated dynamically as we process the signal
avgFilter(lag+1, 1) = mean(y(1:lag+1));
stdFilter(lag+1, 1) = std(y(1:lag+1));

%% Main processing loop
% Process each data point starting after the initialization window
for i = lag+2:length(y)
    
    % Check if current value deviates significantly from expected value
    % "Significantly" means more than (threshold * standard_deviation) away from mean
    if abs(y(i) - avgFilter(i-1)) > threshold * stdFilter(i-1)
        
        % Positive deviation indicates a signal (movement/contraction)
        if (y(i) > avgFilter(i-1))
            signals(i) = 1;  % Mark as signal detected
        else
            % Negative deviation is typically noise, not a true signal
            signals(i) = 0;  % Mark as no signal
        end
        
        % Apply influence factor to the filtered signal
        % This prevents detected peaks from corrupting future thresholds
        % Lower influence = detected signals have less impact on moving statistics
        filteredY(i) = influence * y(i) + (1 - influence) * filteredY(i-1);
        
    else
        % Value is within expected range - no signal detected
        signals(i) = 0;
        
        % Use actual value for filter since it's not an outlier
        filteredY(i) = y(i);
    end
    
    % Update moving statistics using the filtered signal
    % This creates a dynamic threshold that adapts to signal characteristics
    avgFilter(i) = mean(filteredY(i-lag:i));
    stdFilter(i) = std(filteredY(i-lag:i));
end

end
