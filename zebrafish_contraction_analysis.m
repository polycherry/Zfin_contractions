%% ZEBRAFISH EMBRYO CONTRACTION ANALYSIS - MAIN SCRIPT
% =========================================================================
% This is the main script for automated detection and analysis of spontaneous
% contractions in zebrafish embryos from video recordings.
%
% OVERVIEW:
%   Analyzes videos of zebrafish embryos to automatically detect and quantify
%   spontaneous contractions. The pipeline includes:
%   1. Image segmentation to identify embryo positions (ROIs)
%   2. Frame-by-frame motion detection via frame subtraction
%   3. Signal processing to convert noisy motion data into binary signals
%   4. Peak detection to identify individual contraction events
%   5. Data export and visualization
%
% METHODOLOGY:
%   - Videos are recorded at 10 fps
%   - First frame is used as reference for ROI detection via circle detection
%   - Frame subtraction detects movement within each embryo ROI
%   - Robust thresholding algorithm converts signals to binary (moving/not moving)
%   - Contractions are identified as peaks in the binary signal
%
% INPUT:
%   Video files must be named in the format: [group][batch][hours]h.avi
%   Example: "wt1 16h.avi" = wildtype, batch 1, 16 hours post fertilization
%
%   Groups:
%     - wt  = wildtype
%     - mut = mutant
%     - Het = heterozygous
%
% OUTPUT:
%   - freq.xlsx       : Contraction frequency data for each embryo
%   - locmatrix.xlsx  : Timing of each individual contraction event
%   - rawdata.xlsx    : Raw frame subtraction signals
%   - Figure files    : Visualization of detected contractions (.jpg and .fig)
%
% PARAMETERS (Optimized through systematic validation):
%   minpeakdist = 10    : Minimum frames between contractions (1 second)
%   zscore      = 10    : Z-score threshold for signal detection
%   timelag     = 30    : Lag for moving statistics (3 seconds window)
%   smoothing   = 0.5   : Influence factor for threshold adaptation
%
% VALIDATION:
%   Parameters were validated against manual counting of 12 random embryos
%   testing all permutations of plausible parameter values. The chosen
%   combination resulted in sum of squared differences = 18 with minimal
%   false positives and false negatives.
%
% REQUIREMENTS:
%   - MATLAB Image Processing Toolbox (for imfindcircles, im2bw, etc.)
%   - Video files in .avi format
%   - ThresholdingAlgo.m function (must be in same directory or path)
%
% USAGE:
%   1. Place all video files in the same directory as this script
%   2. Ensure video naming follows the required format
%   3. Run this script
%   4. Review generated figures for quality control
%   5. Manually remove dead/falsely detected embryos from output files
%
% NOTES:
%   - Circle detection assumes chorion is nearly perfectly circular
%   - ROI positions remain constant throughout video (no tracking needed)
%   - Dead embryos (visible coagulation) should be manually removed
%   - False detections from Otsu thresholding should be manually filtered
%
% Author: Vasishta Polisetty
% Developed under guidance of: Urvashi Jha (PhD student)
% Lab: Vatsala Thirumalai Lab
% Institution: NCBS (National Centre for Biological Sciences), Bengaluru, India
% Date: 2017
% Updated: February 2026 for GitHub release
% =========================================================================

%% PARAMETER SETUP
% These parameters were optimized through systematic validation
fprintf('Initializing analysis parameters...\n');

minpeakdist = 10;   % Minimum distance between peaks (frames)
                     % = 1 second at 10 fps (documented contraction duration)
                     
zscore = 10;         % Z-score threshold for ThresholdingAlgo
                     % Number of standard deviations for signal detection
                     
timelag = 30;        % Lag for moving average/std calculation (frames)
                     % = 3 second window at 10 fps
                     
smoothing = 0.5;     % Influence/smoothing coefficient (0-1)
                     % Controls how detected signals affect future thresholds

fprintf('Parameters set:\n');
fprintf('  Min peak distance: %d frames (%.1f seconds)\n', minpeakdist, minpeakdist/10);
fprintf('  Z-score threshold: %d\n', zscore);
fprintf('  Time lag: %d frames (%.1f seconds)\n', timelag, timelag/10);
fprintf('  Smoothing: %.1f\n\n', smoothing);

%% INITIALIZE OUTPUT DATA STRUCTURES
% Cell arrays for flexible data storage
freq = cell(0);       % Frequency data: [video_name, embryo_code, embryo_num, num_peaks, frequency]
locmatrix = cell(0);  % Location data: [video_name, embryo_code, embryo_num, peak_location, peak_time]

% Counters for tracking data storage positions
counterA = 1;  % Counter for freq matrix rows
counterB = 1;  % Counter for locmatrix rows
counterC = 1;  % Counter for raw data matrix rows

fprintf('Starting batch processing...\n\n');

%% MAIN PROCESSING LOOPS
% Triple nested loop structure:
%   - Outer loop: Iterate through experimental groups (wt, mut, Het)
%   - Middle loop: Iterate through batches within each group
%   - Inner loop: Iterate through developmental time points (16-27 hours)

for groupnum = 1:3      % Loop through experimental groups
    
    for batch = 1:5     % Loop through batches
        
        % Construct group name string based on group number
        if groupnum == 1
            group = char(strcat('wt', num2str(batch), {' '}));
            groupname = 'wildtype';
        elseif groupnum == 2
            group = char(strcat('mut', num2str(batch), {' '}));
            groupname = 'mutant';
        elseif groupnum == 3
            group = char(strcat('Het', num2str(batch), {' '}));
            groupname = 'heterozygous';
        end
        
        fprintf('Processing %s, batch %d...\n', groupname, batch);
        
        for counter1 = 16:27  % Loop through time points (hours post fertilization)
            
            %% Generate unique code for this embryo video
            % Format: GGBBBHHEE
            % GG = group (01, 02, 03)
            % BBB = batch (00001-00005)
            % HH = hours (16-27)
            % EE = embryo number (added later, 01-99)
            code = groupnum*100000 + batch*10000 + counter1*100;
            
            %% Construct video filename and attempt to load
            str = char(strcat({group}, num2str(counter1), 'h.avi'));
            
            % Try-catch block to handle missing video files gracefully
            try
                vid = VideoReader(str);
            catch err
                fprintf('  Video not found: %s (skipping)\n', str);
                continue  % Skip to next time point if video doesn't exist
            end
            
            % Extract video properties
            frames = vid.NumberOfFrames;
            dur = vid.Duration;  % Duration in seconds
            
            fprintf('  Analyzing %s (%d frames, %.1f sec)...\n', str, frames, dur);
            
            %% IMAGE SEGMENTATION - Detect embryo positions using first frame
            % The first frame serves as reference for ROI detection
            % Embryos remain in approximately the same position throughout
            
            % Import first frame and convert to grayscale
            i = rgb2gray(read(vid, 1));
            [row, col] = size(i);
            
            % Convert to binary using Otsu's method and remove small noise
            % graythresh(i) automatically determines optimal threshold
            % bwareaopen removes connected components smaller than 50 pixels
            bw = bwareaopen(im2bw(i, graythresh(i)), 50);
            
            % Detect circular embryos (chorions) using Hough transform
            % Radius range: 30-45 pixels (typical chorion size)
            % High sensitivity (0.9) to detect most embryos
            % Two-stage method for better accuracy
            [centers, radii] = imfindcircles(bw, [30 45], ...
                                            'Sensitivity', 0.9, ...
                                            'Method', 'twostage');
            
            fprintf('    Detected %d embryos\n', size(centers, 1));
            
            %% VIDEO LOADING - Import all frames into 3D matrix
            % Structure: mov(height, width, frame_number)
            
            % Pre-allocate 3D matrix for efficiency
            mov = zeros(vid.Height, vid.Width, frames, 'uint8');
            
            % Load all frames into memory
            % Only keep first channel if RGB (convert to grayscale)
            counter2 = 1;
            for k = 1:frames
                iim = read(vid, counter2);
                iim = iim(:,:,1);  % Take first channel only
                mov(:,:,counter2) = iim;
                counter2 = counter2 + 1;
            end
            
            %% INITIALIZE PROCESSING MATRICES
            % abcd  : Binary signals (after thresholding) - temporary
            % abcde : Raw frame subtraction signals - temporary
            % abcdef: Permanent storage of all raw data for export
            
            abcd = zeros(size(radii, 1), (ceil(frames/2)+6));
            abcde = zeros(size(radii, 1), (ceil(frames/2)+6));
            
            %% FRAME SUBTRACTION LOOP - Process each detected embryo
            for counter3 = 1:size(radii, 1)
                
                % Create circular mask for current embryo ROI
                % This isolates movement detection to within the chorion
                circleImage = false(row, col);
                [xx, yy] = meshgrid(1:col, 1:row);
                circleImage((xx - centers(counter3,1)).^2 + ...
                           (yy - centers(counter3,2)).^2 <= radii(counter3).^2) = true;
                
                counter4 = 5;  % Starting column for data storage
                
                % Store metadata for this embryo in permanent matrix
                abcdef(counterC, 1) = code + counter3;           % Unique embryo ID
                abcdef(counterC, 2) = dur - (timelag/10);       % Effective duration
                % (subtract timelag since no peaks searched in initial window)
                
                %% Frame-by-frame subtraction within ROI
                for counter5 = 2:1:frames
                    
                    % Subtract consecutive frames to detect motion
                    maskedImageprev = (mov(:,:, counter5)) - (mov(:,:, counter5-1));
                    
                    % Apply circular mask - zero out everything outside embryo
                    maskedImageprev(~circleImage) = 0;
                    
                    % Calculate mean intensity of difference and square it
                    % Squaring: 1) Makes all values positive
                    %          2) Amplifies strong signals over noise
                    signal_value = mean2(maskedImageprev).^2;
                    
                    % Store in both temporary and permanent matrices
                    abcde(counter3, counter4) = signal_value;
                    abcdef(counterC, counter4) = signal_value;
                    
                    counter4 = counter4 + 1;
                end
                
                counterC = counterC + 1;  % Move to next row for next embryo
            end
            
            %% SIGNAL PROCESSING - Convert raw signals to binary
            fprintf('    Processing signals with thresholding algorithm...\n');
            
            % Apply robust thresholding algorithm to each embryo's signal
            for xc = 1:size(radii, 1)
                [abcd(xc,:), ~, ~] = ThresholdingAlgo(abcde(xc,:), ...
                                                      timelag, ...
                                                      zscore, ...
                                                      smoothing);
            end
            
            %% PEAK DETECTION - Identify individual contraction events
            fprintf('    Detecting contraction peaks...\n');
            
            for counter6 = 1:size(radii, 1)
                
                % Find peaks in binary signal
                % MinPeakDistance prevents detecting same contraction multiple times
                [peaks, locs] = findpeaks(abcd(counter6, :), ...
                                         'MinPeakDistance', minpeakdist);
                
                % Store frequency data in output matrix
                freq{counterA, 1} = char(strcat(group, num2str(counter1)));  % Video name
                freq{counterA, 2} = code + counter6;                          % Unique ID
                freq{counterA, 3} = counter6;                                 % Embryo number
                freq{counterA, 4} = length(peaks);                            % Total contractions
                freq{counterA, 5} = length(peaks) / (dur - (timelag/10));   % Frequency (Hz)
                
                counterA = counterA + 1;
                
                % Store individual peak locations/times
                if isempty(locs)
                    % No peaks detected - store zero entry
                    locmatrix{counterB, 1} = char(strcat(group, num2str(counter1)));
                    locmatrix{counterB, 2} = code + counter6;
                    locmatrix{counterB, 3} = counter6;
                    locmatrix{counterB, 4} = 0;      % Peak location (frames)
                    locmatrix{counterB, 5} = 0;      % Peak time (seconds)
                    
                    counterB = counterB + 1;
                else
                    % Store each peak location separately
                    for counter7 = 1:length(locs)
                        locmatrix{counterB, 1} = char(strcat(group, num2str(counter1)));
                        locmatrix{counterB, 2} = code + counter6;
                        locmatrix{counterB, 3} = counter6;
                        locmatrix{counterB, 4} = locs(counter7);           % Frame number
                        locmatrix{counterB, 5} = locs(counter7) / 10;      % Time (seconds)
                        
                        counterB = counterB + 1;
                    end
                end
            end
            
            % Status update
            status = sprintf('    Completed %d - detected %d embryos with %d total contractions', ...
                           code + counter1, size(radii, 1), sum(cell2mat(freq(end-size(radii,1)+1:end, 4))));
            fprintf('%s\n', status);
            
            %% VISUALIZATION - Create diagnostic figure
            h1 = figure(1);
            set(h1, 'Position', [100, 100, 1200, 800]);  % Large figure for clarity
            
            counter8 = 1;
            
            % Create subplot for each embryo
            for counter9 = 1:length(radii)
                
                emb = char(strcat('e', num2str(counter9)));  % Embryo label
                
                % Find peaks for this embryo
                [pks1, locs1] = findpeaks(abcd(counter9, :), ...
                                         'MinPeakDistance', minpeakdist);
                
                %% Plot raw signal with detected peaks
                if (counter9 < length(radii))
                    % Not last plot - omit x-axis labels to save space
                    subplot(length(radii), 2, counter8);
                    
                    % Plot signal
                    plot(abcde(counter9,:));
                    hold on;
                    
                    % Mark detected peaks
                    if (length(locs1) ~= 0)
                        plot(locs1, 0, 'r*', 'markerfacecolor', [1 0 0]);
                    end
                    
                    % Format plot
                    set(gca, 'Xtick', 304:300:1804);
                    set(gca, 'Xticklabel', '');  % No labels except last plot
                    set(gca, 'Ytick', 0);
                    set(gca, 'Yticklabel', emb);  % Label with embryo number
                    grid on;
                    hold off;
                    
                else
                    % Last plot - include x-axis labels (time in minutes)
                    subplot(length(radii), 2, counter8);
                    
                    plot(abcde(counter9,:));
                    hold on;
                    
                    if (length(locs1) ~= 0)
                        plot(locs1, 0, 'r*', 'markerfacecolor', [1 0 0]);
                    end
                    
                    % Format with time labels
                    set(gca, 'Xtick', 304:300:1804);
                    set(gca, 'Xticklabel', {'.5', '1', '1.5', '2', '2.5', '3min'});
                    set(gca, 'Ytick', 0);
                    set(gca, 'Yticklabel', emb);
                    grid on;
                    hold off;
                end
                
                counter8 = counter8 + 2;
            end
            
            %% Display reference image with labeled embryos
            % Right column shows first frame with detected circles
            subplot(length(radii), 2, [2:2:((length(radii))*2)]);
            
            % Create labels for each embryo
            names = cell(length(radii), 1);
            for run = 1:length(radii)
                names{run} = strcat('e', num2str(run));
            end
            
            % Annotate image with embryo numbers
            i = insertText(i, centers, names);
            imshow(i)
            
            % Draw circles around detected embryos
            h = viscircles(centers, radii);
            title(sprintf('%s - Otsu threshold = %.3f', vid.Name, graythresh(i)));
            
            %% Save figure for manual quality control
            saveas(h1, char(strcat({group}, num2str(counter1))), 'jpg');
            saveas(h1, char(strcat({group}, num2str(counter1))), 'fig');
            
            close(h1);  % Close figure to free memory
            
        end  % End time point loop
        
    end  % End batch loop
    
    fprintf('\n');
    
end  % End group loop

%% EXPORT RESULTS
fprintf('\nExporting results...\n');

% Export frequency data
fprintf('  Writing freq.xlsx...\n');
xlswrite('freq.xlsx', freq, 1);

% Export peak location/timing data
fprintf('  Writing locmatrix.xlsx...\n');
xlswrite('locmatrix.xlsx', locmatrix, 1);

% Export raw frame subtraction signals
fprintf('  Writing rawdata.xlsx...\n');
xlswrite('rawdata.xlsx', abcdef, 1);

fprintf('\n=== Analysis Complete ===\n');
fprintf('Total embryos processed: %d\n', counterA-1);
fprintf('Total contractions detected: %d\n', counterB-1);
fprintf('\nPlease review generated figures for quality control.\n');
fprintf('Manually remove dead or falsely detected embryos from output files.\n');
