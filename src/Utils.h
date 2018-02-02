#ifndef UTILS
#define UTILS

#include <vector>
#include <iostream>
#include <cmath>
#include "Vector.h"

void generateData(std::vector<Vector>& x, std::vector<double>& y, int s);

void generateNoisyData(std::vector<Vector>& x, std::vector<double>& y, int s);

void generateData2(std::vector<Vector>& x, std::vector<double>& y, int s);

void generateData3(std::vector<Vector>& x,std::vector<double>& y, int s);

#endif
