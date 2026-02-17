# Zebrafish Embryo Contraction Analysis

Automated detection and quantification of spontaneous contractions in zebrafish embryos from video recordings using MATLAB.

## Overview

This project provides an automated pipeline for analyzing spontaneous contractions in zebrafish embryos during early development. Traditional methods rely on manual counting or proprietary software, limiting sample size and temporal precision. This automated approach enables:

- **Higher throughput**: Process multiple embryos automatically
- **Improved precision**: Frame-by-frame analysis at 10 fps
- **Rich time-series data**: Exact timing of each contraction, not just frequency
- **Reproducible analysis**: Validated algorithm eliminates observer bias

## Quick Start

### Prerequisites

- MATLAB (R2016b or later)
- Image Processing Toolbox
- Video files in `.avi` format

### Installation

```bash
git clone https://github.com/vasishtapolisetty/Zfin_contractions.git
cd Zfin_contractions
```

### Usage

1. **Name your videos**: `[group][batch] [hours]h.avi`
   - Example: `wt1 16h.avi` = wildtype, batch 1, 16 hours post-fertilization

2. **Place videos** in the same directory as the MATLAB scripts

3. **Run the analysis**:
   ```matlab
   zebrafish_contraction_analysis
   ```

4. **Review outputs**:
   - `freq.xlsx` - Contraction frequency for each embryo
   - `locmatrix.xlsx` - Timing of individual contractions
   - `rawdata.xlsx` - Raw frame subtraction signals
   - Diagnostic figures (.jpg/.fig) for quality control

5. **Quality control**: Manually review figures and remove dead/falsely detected embryos

## How It Works

### 1. Image Segmentation
- First frame converted to binary using Otsu thresholding
- Hough transform circle detection identifies embryo positions
- ROIs remain constant throughout video (chorions are stationary)

### 2. Video Processing
- Frame-by-frame subtraction detects motion within each ROI
- Squared differences amplify signal over noise
- Each embryo processed independently

### 3. Signal Processing (ThresholdingAlgo)
The core innovation uses dynamic z-score thresholding:
- Constructs moving mean and standard deviation
- Uses relative thresholds (z-scores), not absolute values
- Prevents detected signals from corrupting future thresholds
- Converts continuous signal to binary (moving/not moving)

### 4. Peak Detection
- `findpeaks()` identifies contraction events from binary signal
- Minimum peak distance = 1 second (documented contraction duration)

## Algorithm Parameters

Optimized through testing 1,200 parameter combinations against manual counting of 12 embryos:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `timelag` | 30 frames | Moving statistics window (3 seconds) |
| `zscore` | 10 | Z-score threshold for detection |
| `smoothing` | 0.5 | Influence factor (0-1) |
| `minpeakdist` | 10 frames | Min separation between peaks (1 second) |

### Validation Results

Parameters were validated against manual counting of 12 random embryos:

| Embryo ID | Manual Count | Automated Count | False Positives | False Negatives |
|-----------|--------------|-----------------|-----------------|-----------------|
| 212604    | 19           | 19              | 0               | 0               |
| 212625    | 21           | 20              | 0               | 1               |
| 322303    | 13           | 12              | 1               | 2               |
| 322601    | 18           | 17              | 0               | 1               |
| 321901    | 0            | 0               | 0               | 0               |
| 321902    | 0            | 0               | 0               | 0               |
| 322602    | 18           | 18              | 0               | 0               |
| 322605    | 14           | 14              | 0               | 0               |
| 122605    | 14           | 14              | 0               | 0               |
| 111909    | 21           | 21              | 0               | 0               |
| 312002    | 4            | 4               | 0               | 0               |
| 111915    | 32           | 28              | 0               | 5               |
| **Total** | **174**      | **167**         | **1**           | **9**           |

**Overall Performance:**
- **Accuracy**: 96.0% (167/174 correct detections)
- **Precision**: 99.4% (166/167 detected were true positives)
- **Sensitivity**: 94.8% (165/174 actual contractions detected)
- **Sum of squared differences**: 18 (best among 1,200 parameter combinations tested)

### Adjusting Parameters

Edit lines 15-18 in `zebrafish_contraction_analysis.m`:

```matlab
minpeakdist = 10;   # Min frames between peaks
zscore = 10;        # Z-score threshold
timelag = 30;       # Moving average window
smoothing = 0.5;    # Influence factor
```

**Common adjustments:**
- **Too many false positives**: Increase `zscore` to 15
- **Missing weak contractions**: Decrease `zscore` to 7
- **Different frame rate**: Scale `minpeakdist` and `timelag` proportionally

## Output Format

### freq.xlsx
| Video Name | Embryo Code | Embryo # | Total Contractions | Frequency (Hz) |
|------------|-------------|----------|-------------------|----------------|
| wt1 16 | 111601 | 1 | 23 | 0.128 |

### locmatrix.xlsx
| Video Name | Embryo Code | Embryo # | Frame | Time (sec) |
|------------|-------------|----------|-------|------------|
| wt1 16 | 111601 | 1 | 150 | 15.0 |
| wt1 16 | 111601 | 1 | 225 | 22.5 |

**Embryo code format**: `GGBBBHHEE`
- GG = Group (01=wt, 02=mut, 03=Het)
- BBB = Batch (1-5)
- HH = Hours (16-27)
- EE = Embryo number

## Troubleshooting

### No embryos detected
- Adjust circle detection radius in `imfindcircles`: default is `[30 45]`
- Increase sensitivity: try `0.95`
- Check that first frame has good contrast

### Memory errors
- Process fewer videos at once
- Clear workspace between batches: `clear mov abcd abcde`

### Wrong frame rate
Calculate new parameters:
```matlab
fps = 15;  % Your frame rate
minpeakdist = round(1.0 * fps);   % 1 second
timelag = round(3.0 * fps);        % 3 seconds
```

## Technical Details

### Mathematical Foundation

**Z-score calculation:**
```
z = (y - μ) / σ
```
Signal detected if `|z| > zscore` where μ and σ are moving statistics.

**Moving statistics:**
```
μ(t) = mean(y[t-lag:t])
σ(t) = std(y[t-lag:t])
```

**Influence factor:**
```
y*(t) = smoothing·y(t) + (1-smoothing)·y*(t-1)  [if signal detected]
y*(t) = y(t)                                      [otherwise]
```

### Video Requirements

**Recommended settings:**
- Frame rate: 10 fps
- Duration: 3-5 minutes
- Resolution: 640×480 or higher
- Lighting: Even, transmitted light
- Format: .avi with minimal compression

**Recording tips:**
- Focus on chorion edges for best circle detection
- Avoid overlapping embryos
- Keep chorions fully in frame
- Minimize bubbles and debris

## Citation

If you use this code, please cite:

```
Polisetty, V. (2017). Automated Detection of Zebrafish Embryo Contractions.
Developed at Vatsala Thirumalai Lab, NCBS Bengaluru.
GitHub: https://github.com/vasishtapolisetty/Zfin_contractions
```

## Acknowledgments

Developed in **2017** at the **National Centre for Biological Sciences (NCBS), Bengaluru, India**.

**Special thanks to:**
- **Urvashi Jha** (PhD student) - Mentorship and validation methodology
- **Vatsala Thirumalai Lab** - Resources and research environment

## Contributing

Contributions and suggestions are welcome! Please open an issue to discuss proposed changes.

## Contact

- **GitHub Issues**: [Report bugs or request features](https://github.com/polycherry/Zfin_contractions/issues)
- **Email**: vasishtapolisetty@gmail.com

## License

MIT License - see [LICENSE](LICENSE) file for details.

---

**Original Development**: 2017, NCBS Bengaluru  
**GitHub Release**: February 2026
