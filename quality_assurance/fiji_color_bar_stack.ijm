// Parameters
barWidth = 20;            // Colorbar width in pixels
barHeightFactor = 0.9;    // Fraction of image height for colorbar height
lutName = "Fire";         // Use Fiji's 'Fire' LUT for colorbar

// Get current image
id = getImageID();

stackSize = nSlices;
width = getWidth();
height = getHeight();
barHeight = round(height * barHeightFactor);

// Get the actual min and max values
getMinAndMax(min, max);

// Calculate new max for colorbar = min + 0.75*(max-min)
barMax = min + 0.5 * (max - min);

// Create colorbar image (float)
newImage("Colorbar", "32-bit black", barWidth, barHeight, 1);
colorbarID = getImageID();

// Fill colorbar (bottom = min, top = barMax)
for (y = 0; y < barHeight; y++) {
    val = min + (barMax - min) * y / (barHeight - 1);
    for (x = 0; x < barWidth; x++) {
        setPixel(x, barHeight - 1 - y, val);
    }
}

// Apply LUT (display only)
selectImage(colorbarID);
run(lutName);

// Create new stack with the colorbar at right, in float format
newWidth = width + barWidth;
newStackTitle = getTitle() + "-with-colorbar";
newImage(newStackTitle, "32-bit black", newWidth, height, stackSize);

for (i = 1; i <= stackSize; i++) {
    // Paste original slice
    selectImage(id);
    setSlice(i);
    run("Copy");
    selectImage(newStackTitle);
    setSlice(i);
    makeRectangle(0, 0, width, height);
    run("Paste");

    // Paste colorbar
    selectImage(colorbarID);
    makeRectangle(0, 0, barWidth, barHeight);
    run("Copy");
    selectImage(newStackTitle);
    setSlice(i);
    yOffset = round((height - barHeight) / 2);
    makeRectangle(width, yOffset, barWidth, barHeight);
    run("Paste");
}

// Set display range and LUT for the new stack
selectImage(newStackTitle);
setMinAndMax(min, barMax);
run(lutName);

// Cleanup
selectImage(colorbarID); close();