#include "Utils.h"

void generateData(std::vector<Vector>& x,std::vector<double>& y, int s) {
    srand(time(NULL));
    int lower_bound_x = -4, upper_bound_x = 4;
    int lower_bound_y = -1, upper_bound_y = 1;
    for (int i = 0; i<s; i++) {
        double input = ((double)rand() / (double)RAND_MAX) * 6.28 - 3.14;
        Vector xi(1);
        xi(0) = (input - lower_bound_x) / (upper_bound_x - lower_bound_x);
        x.push_back(xi);
        y.push_back((std::sin(input) - lower_bound_y) / (upper_bound_y - lower_bound_y));
    }
}
