# Kymograph fibril tip tracking

Tools for extracting the position of a growing fibril tip from kymographs and testing whether its growth rate changes over time.

The workflow has two steps:

1. **Fiji/ImageJ macro** (`export_kymograph_matrix_multiple.ijm`) exports kymograph images as plain-text pixel matrices.
2. **R script** (`kymograph_edge_extraction_annotated.R`) reads a matrix, detects the tip position in every frame, and runs change-point analysis on the resulting trace.

## Repository contents

```
.
├── export_kymograph_matrix_multiple.ijm     # Fiji macro: kymographs -> text matrices
├── kymograph_edge_extraction_annotated.R    # R script: tip tracking + change-point analysis
├── data/                                    # exported matrices (e.g. fibril15_matrix.txt)
├── renv.lock                                # R package versions
├── renv/                                    # renv set-up files (activate.R, settings.json)
├── .Rprofile                                # activates renv when the project is opened
└── README.md
```

The `data/` folder is where the R script looks for matrices by default. Adjust the folder layout to match your repository.

## Requirements

- [Fiji](https://fiji.sc/) (ImageJ) for the export step
- R, with the packages recorded in `renv.lock` (the lockfile was created with R 4.2.2)
  - `tidyverse`
  - `changepoint`
  - `segmented`

To restore the R environment, open the project in RStudio (or start R in the project folder) and run:

```r
renv::restore()
```

## Step 1: export kymographs from Fiji

The macro starts from kymograph images that are already open in Fiji. It does not generate kymographs.

Kymograph orientation must be:

- **vertical axis = time** (one row per frame)
- **horizontal axis = position along the fibril** (pixels)

To export:

1. Open all the kymograph images you want to export in Fiji.
2. Run `export_kymograph_matrix_multiple.ijm` (Plugins > Macros > Run...).
3. Choose an output folder when prompted.

Each open image is saved as a tab-delimited text file (Fiji's "Text Image" format) named after its window title, with the extension replaced by `_matrix.txt`. For example, `fibril_03.tif` becomes `fibril_03_matrix.txt`. The images are left open afterwards.

Each exported file has no header, rows are frames (first row = first frame), columns are pixel positions, and brighter values mean fibril signal.

## Step 2: analyse a matrix in R

The script processes **one fibril (one matrix file) per run**. Set the parameters below, then run the script from top to bottom.

### Parameters to set

| Parameter | Where | Meaning |
|---|---|---|
| `matrix_file` | Section 1 | Path to the exported matrix, e.g. `data/fibril15_matrix.txt` |
| `n_background_frames` | Section 2 | Number of frames at the start of the movie that contain background only (before any fibril tip is visible) |
| `anchor_columns` | Section 2 | Position columns next to the fixed anchor/nucleation point that are covered by fibril for the rest of the movie once growth has started |
| `background_side` | Section 2 | `"right"` if fibrils grow towards higher columns, `"left"` if they grow towards lower columns |
| `px_size_nm` | Section 3 | Nanometres per pixel along the position axis |
| `time_per_frame` | Section 3 | Minutes per frame |

The values in the script (40 frames, columns 1:5, 160 nm per pixel, 2 min per frame) are example values and need to be set for your own acquisition.

### What the script does

1. **Load** the matrix.
2. **Set a threshold.** The background level is the median intensity of the initial background-only frames. The fibril level is the median intensity of the anchor columns in the frames after growth starts. The threshold is the midpoint between the two, calculated once and applied to every frame.
3. **Detect the tip** in each frame by finding where the intensity profile crosses the threshold, with sub-pixel precision (linear interpolation). Frames where the tip cannot be located unambiguously are dropped.
4. **Build the trace:** tip position (nm) against time (min), and plot it.
5. **Change-point analysis:**
   - `changepoint::cpt.mean` finds jumps in tip position (exploratory).
   - `segmented::selgmented` fits a piecewise-linear model with 0 to 3 breakpoints, choosing the number by BIC, to detect changes in growth rate.
6. **Export** an SVG of the tip trace with a linear fit.
7. **Report** breakpoint locations with confidence intervals and the growth rate (slope, nm/min) of each segment.

### Outputs

- Plots in the plot pane: tip trace, `cpt.mean` change-points, segmented fit
- An SVG saved in the working directory (the filename is set in section 6 of the script)
- Console output: background, fibril and threshold intensities, the segmented model summary, breakpoint confidence intervals and segment slopes

If no breakpoint is selected, the script reports the slope of a single linear fit instead.

### Checking the output

- After the threshold is calculated, check that the global contrast printed to the console is clearly positive. If it is small or negative, check `n_background_frames` and `anchor_columns`.
- Inspect the tip trace plot by eye before interpreting the change-point results.

## Notes and limitations

- Fibrils are assumed to be anchored at a fixed point, with only the tip moving.
- Tip detection uses the last pixel above threshold, so an isolated bright noise pixel beyond the real tip can shift the detected edge.
- The first detected point is excluded before the segmented fit, because it is often not detected correctly.
- `cpt.mean` uses a Normal cost that assumes unit variance. A trace in nm that has not been scaled may give many change-points, so treat this analysis as exploratory.
