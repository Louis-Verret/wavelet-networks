#ifndef UTILS
#define UTILS

#include <vector>
#include <iostream>
#include <cmath>
#include "Vector.h"

/* Generate sinus data for training and evaluation */
void generateData(std::vector<Vector>& x, std::vector<double>& y, int s);

/* Generate noisy sinus data for training and evaluation */
void generateNoisyData(std::vector<Vector>& x, std::vector<double>& y, int s);

/* Generate complex 1D data as described in the paper, for training and evaluation */
void generateData2(std::vector<Vector>& x, std::vector<double>& y, int s);

/* Generate 2D data as described in the paper, for training and evaluation */
void generateData3(std::vector<Vector>& x,std::vector<double>& y, int s);

#endif
