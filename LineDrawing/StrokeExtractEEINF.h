/*
Based on the paper
Stroke Extraction in Cartoon Images Using Edge-Enhanced Isotropic Nonlinear Filter
VRCAI '10

[Assuming the image has already been smoothed]
*/

int StrokeExtract(unsigned char * image,
    int height, int width, int channel,
    int radius, double epsilon);