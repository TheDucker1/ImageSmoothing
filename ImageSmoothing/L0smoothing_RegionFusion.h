/*
The Code is created based on the method described in the following paper
"Fast and Effective L0 Gradient Minimization by Region Fusion", 
Rang M. H. Nguyen, Michael S. Brown
ICCV 2015
*/

int L0smoothing(unsigned char * image, 
    int height, int width, int channel,
    unsigned char * mask, float lambda);