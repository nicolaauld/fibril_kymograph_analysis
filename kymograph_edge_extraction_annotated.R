# kymograph_edge_extraction.R
#
# PURPOSE
#   Tracks the growing tip of a single fibril in a kymograph and tests
#   whether its growth rate changes over time.
#
# WORKFLOW
#   1. Load a kymograph exported from Fiji as a plain-text pixel matrix.
#   2. Estimate background and fibril intensity levels, and set a single
#      threshold halfway between them.
#   3. For each frame, find where the intensity profile crosses that
#      threshold (the fibril tip), with sub-pixel precision.
#   4. Convert the tip position to nm and time to minutes, and plot it.
#   5. Run change-point analysis on the trace:
#        A. changepoint::cpt.mean   - jumps in tip position (exploratory)
#        B. segmented::selgmented   - changes in growth rate (slope)
#   6. Export a figure and report breakpoints and segment-wise growth rates.
#
# INPUT FILE REQUIRED
#   One tab-delimited .txt file per fibril, exported from Fiji
#   (File > Save As > Text Image, or via export_kymograph_matrix.ijm):
#     - no header row, no row/column names
#     - rows    = frames (time), first row = first frame
#     - columns = pixels along the fibril (position)
#     - bright pixels = fibril present, dark pixels = background
#   The fibril is assumed to be anchored at a fixed point (the columns
#   nearest the anchor are always covered by fibril once growth starts)
#   and to grow towards the right (higher column numbers). If your
#   fibrils grow towards the left, set background_side = "left" in
#   section 2.
#
# R PACKAGES REQUIRED
#   tidyverse, changepoint, segmented
#   (if using renv: renv::restore() installs the versions in renv.lock)
#
# OUTPUT
#   - Plots in the RStudio plot pane (tip trace, change-points, segmented fit)
#   - An SVG of the tip trace with a linear fit, saved in the working directory
#   - Console output: reference intensity levels, segmented model summary,
#     breakpoint confidence intervals and segment slopes
#
# PARAMETERS YOU MUST SET (marked "EDIT")
#   matrix_file, n_background_frames, anchor_columns,
#   px_size_nm, time_per_frame

library(tidyverse)
library(changepoint)
library(segmented)

# ---- 1. Load the exported matrix -------------------------------------------

# EDIT: path to the exported kymograph matrix 
matrix_file <- "data/fibril15_matrix.txt"

# header = FALSE because the Fiji text export has no column names.
mat <- as.matrix(
  read.table(matrix_file, header = FALSE, sep = "\t"))

n_frames    <- nrow(mat)   # number of time points
n_positions <- ncol(mat)   # number of pixels along the fibril axis

# ---- 2. Per-frame edge detection -------------------------------------------

# EDIT these two values for your movie:
# - n_background_frames: number of frames at the START of the movie that are
#   reliably background-only (before any fibril tip is visible).
# - anchor_columns: position columns that are reliably covered by fibril for
#   the rest of the movie once growth has started. These are the columns next
#   to the fixed nucleation/anchor point, which does not move (only the tip
#   does). Check them against your kymograph image.
n_background_frames <- 40    # set to a frame count you trust
anchor_columns      <- 1:5   # columns nearest the anchor point

# Reference intensity levels, each computed ONCE for the whole kymograph:
#   - background: median of all pixels in the initial background-only frames
#     (no fibril exists yet, so these pixels are background by definition).
#   - fibril: median of the anchor columns in the frames after growth starts
#     (the base of the fibril stays put, so these pixels are always fibril).
# Using fixed levels means every frame is judged against the same threshold,
# so noise in an individual frame cannot shift the threshold itself.
bg_level <- median(mat[1:n_background_frames, ], na.rm = TRUE)

fg_level <- median(
  mat[(n_background_frames + 1):n_frames, anchor_columns],
  na.rm = TRUE
)

# Threshold = halfway between background and fibril intensity.
threshold       <- (bg_level + fg_level) / 2
contrast_global <- fg_level - bg_level

cat("bg_level:", bg_level, "\n")
cat("fg_level:", fg_level, "\n")
cat("threshold:", threshold, "\n")
cat("global contrast:", contrast_global, "\n")
# Sanity check: contrast_global should be clearly positive and reasonably
# large. If it is small or negative, check n_background_frames and
# anchor_columns (e.g. anchor_columns overlapping background, or growth
# already starting within the first n_background_frames).

# find_edge(): locates the fibril tip in one frame (one row of the matrix).
# It finds the last (or first) pixel above the threshold, then linearly
# interpolates between that pixel and its neighbour to estimate where the
# profile crosses the threshold, giving a sub-pixel position.
#
# Arguments:
#   profile         - pixel intensities for one frame
#   threshold       - the fixed threshold calculated above
#   background_side - which end of the row is dark background:
#                     "right" if the fibril grows towards higher columns,
#                     "left" if it grows towards lower columns
#
# Returns NA when the tip cannot be located unambiguously: the whole row is
# below threshold (no fibril yet), or the profile is still above threshold at
# the far edge of the image (the tip has left the field of view or never
# crosses). These frames are dropped in section 3.
#
# Note: because it uses the last pixel above threshold, an isolated bright
# noise pixel beyond the real tip can shift the detected edge.
find_edge <- function(profile, threshold, background_side = c("right", "left")) {
  background_side <- match.arg(background_side)

  above <- profile > threshold

  if (background_side == "right") {
    # fibril (bright) on the left, background (dark) on the right;
    # find the last position still above threshold
    if (!any(above)) return(NA_real_)             # whole row below threshold
    idx <- max(which(above))
    if (idx == length(profile)) return(NA_real_)  # never crosses -> ambiguous frame
  } else {
    if (!any(above)) return(NA_real_)
    idx <- min(which(above))
    if (idx == 1) return(NA_real_)
  }

  # The two pixels either side of the crossing
  if (background_side == "right") {
    x1 <- idx; x2 <- idx + 1
  } else {
    x1 <- idx - 1; x2 <- idx
  }
  y1 <- profile[x1]; y2 <- profile[x2]

  if (y1 == y2) return((x1 + x2) / 2)  # flat segment: avoid division by zero

  # Position where the profile equals the threshold, interpolated between x1 and x2
  x1 + (threshold - y1) * (x2 - x1) / (y2 - y1)
}

# Apply to every frame using the single fixed threshold.
# Change background_side to "left" if your fibrils grow in the opposite direction.
edge_positions <- sapply(seq_len(n_frames), function(i) {
  find_edge(mat[i, ], threshold = threshold, background_side = "right")
})

# ---- 3. Build the position/time tibble -------------------------------------

# EDIT these conversion factors for your acquisition settings:
px_size_nm     <- 160   # nm per pixel along the position axis
time_per_frame <- 2     # minutes per frame (if you use seconds, change the
                        # axis labels and units below to match)

tip_trace <- tibble(
  frame       = seq_len(n_frames),
  time_min    = (frame - 1) * time_per_frame,   # frame 1 = time 0
  position_px = edge_positions,
  position_nm = edge_positions * px_size_nm
) %>%
  filter(!is.na(position_nm))   # drop frames where the tip was not detected
                                # (e.g. before growth starts); time_min keeps
                                # its original values, so the time axis is unaffected

# Plot the extracted trace and check it by eye before any further analysis.
ggplot(tip_trace, aes(time_min, position_nm)) +
  geom_point() +
  geom_line() +
  labs(x = "Time (min)", y = "Tip position (nm)",
       title = "Extracted fibril tip position over time")

# ---- 4. Change-point analysis: jumps in tip position -----------------------

# cpt.mean() finds points where the MEAN of the trace shifts (piecewise-
# constant model), using the PELT algorithm. This is useful for discrete
# jumps in position. It does not model a steadily increasing trace, so it
# is exploratory here; growth-rate changes are tested in section 5.
# The default (Normal) cost assumes unit variance, so if the trace is not
# standardised (e.g. position in nm) it may report many change-points;
# consider scaling the data first.
cpt_mean <- cpt.mean(tip_trace$position_nm, method = "PELT")
plot(cpt_mean, main = "Change-points in tip position (mean)")

# ---- 5. Change-point analysis: changes in growth rate ----------------------

# For a growing fibril the quantity of interest is usually the growth RATE.
# A piecewise-linear (segmented) regression detects points where the slope
# of position vs time changes, i.e. transitions between growth phases.

# Drop the first detected point before fitting. This is often because it is
# not detecting the first one correctly
tip_trace_edit <- tip_trace[-1, ]

# Single straight-line fit: the starting model for selgmented().
lin_model <- lm(position_nm ~ time_min, data = tip_trace_edit)

# selgmented() tests for 0 up to Kmax breakpoints and selects the number
# that gives the best BIC, so the number of breakpoints is not fixed in
# advance. seg.Z names the variable the breakpoints are placed along.
seg_model <- selgmented(lin_model, Kmax = 3, type = "bic", seg.Z = ~time_min)

summary(seg_model)

# Plot the data with the fitted model.
# If 0 breakpoints fit best, selgmented() returns an ordinary lm object
# rather than a "segmented" object, so the two cases are plotted differently.
plot(tip_trace_edit$time_min, tip_trace_edit$position_nm, pch = 16,
     xlab = "Time (min)", ylab = "Tip position (nm)")

if (inherits(seg_model, "segmented")) {
  plot(seg_model, add = TRUE, col = "red", lwd = 2)
} else {
  abline(lin_model, col = "red", lwd = 2)
}

# ---- 6. Figure export ------------------------------------------------------

# Tip position with a single linear fit, saved as an SVG in the working
# directory. EDIT the filename to suit.
linear_graph <- ggplot(tip_trace_edit, aes(x = time_min, y = position_nm)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()

svg(filename = "no_plateau_example_fibril_2.svg", width = 6, height = 5)
print(linear_graph)
dev.off()

# ---- 7. Breakpoints and growth rates ---------------------------------------

# These only apply if at least one breakpoint was found (seg_model is a
# "segmented" object).
if (inherits(seg_model, "segmented")) {
  # Breakpoint locations (in time_min) with confidence intervals
  print(confint(seg_model))

  # Growth rate (slope, nm per min) of each segment
  print(slope(seg_model))
} else {
  cat("No breakpoints selected: growth is described by a single linear fit.\n")
  print(coef(lin_model))   # slope = growth rate in nm per min
}
