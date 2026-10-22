// export_kymograph_matrix.ijm
//
// Exports every currently OPEN kymograph image as a tab-delimited text
// matrix (one .txt per kymograph). Only touches images already open in
// Fiji - does not scan or read anything else from disk.
// Rows = time (frames), Columns = position along fibril (pixels).
// This matches: vertical (long) axis = time, horizontal (narrow) axis = position.
//
// Usage:
//   1. Open all the kymograph images you want to export (e.g. 10 files).
//   2. Run this macro (Plugins > Macros > Run...).
//   3. Choose an output folder when prompted - all matrices are saved there.
//
// Each output file is named after its source image's window title, e.g.
// "fibril_03.tif" -> "fibril_03_matrix.txt", so you can match them back up
// to the right fibril in R. Windows are left open afterward (not closed).

if (nImages == 0) {
    exit("No images open. Please open your kymograph(s) first.");
}

outputDir = getDirectory("Choose a folder to save the matrix files");

n_total = nImages;

for (i = 1; i <= n_total; i++) {
    selectImage(i);
    title = getTitle();
    defaultName = replace(title, "\\.[a-zA-Z0-9]+$", "") + "_matrix.txt";
    savePath = outputDir + defaultName;

    // Save As > Text Image writes the raw pixel matrix, tab-delimited,
    // with rows corresponding to image rows (y) and columns to image columns (x).
    saveAs("Text Image", savePath);
}

print("Done.");
print("Processed: " + n_total + " open image(s)");
print("Output folder: " + outputDir);
print("Rows = frames (time), Columns = position (pixels) in each output file");
