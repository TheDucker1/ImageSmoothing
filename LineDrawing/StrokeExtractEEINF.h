/*
Based on the paper
Stroke Extraction in Cartoon Images Using Edge-Enhanced Isotropic Nonlinear Filter
VRCAI '10

[Using guided filter instead of  bilateral filter as described in the paper]
*/

int StrokeExtract(unsigned char * image,
    int height, int width, int channel,
    double epsilon);