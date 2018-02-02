#include "Utils.h"
#include <random>

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

void generateNoisyData(std::vector<Vector>& x,std::vector<double>& y, int s) {
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0,0.2);
    int lower_bound_x = -4, upper_bound_x = 4;
    int lower_bound_y = -1, upper_bound_y = 1;
    for (int i = 0; i<s; i++) {
        double input = ((double)rand() / (double)RAND_MAX) * 6.28 - 3.14;
        Vector xi(1);
        xi(0) = (input - lower_bound_x) / (upper_bound_x - lower_bound_x);
        x.push_back(xi);
        double noise = distribution(generator);
        y.push_back((std::sin(input) + noise - lower_bound_y) / (upper_bound_y - lower_bound_y));
    }
}

void generateData2(std::vector<Vector>& x,std::vector<double>& y, int s) {
    srand(time(NULL));
    int lower_bound_x1 = -10, upper_bound_x1 = 10;
    int lower_bound_x2 = -10, upper_bound_x2 = 10;
    int lower_bound_y = -95.9, upper_bound_y = 95.9;
    for (int i = 0; i<s; i++) {
        double x1 = ((double)rand() / (double)RAND_MAX) * 20 - 10;
        double x2 = ((double)rand() / (double)RAND_MAX) * 20 - 10;
        Vector xi(2);
        xi(0) = (x1 - lower_bound_x1) / (upper_bound_x1 - lower_bound_x1);
        xi(1) = (x2 - lower_bound_x2) / (upper_bound_x2 - lower_bound_x2);
        x.push_back(xi);
        y.push_back(((x1 * x1 - x2 * x2) * std::sin(0.5 * x1) - lower_bound_y) / (upper_bound_y - lower_bound_y));
    }
}

void generateData3(std::vector<Vector>& x_vec,std::vector<double>& y, int s) {
    srand(time(NULL));
    int lower_bound_x1 = -10, upper_bound_x1 = 10;
    int lower_bound_y = -10, upper_bound_y = 10;
    for (int i = 0; i<s; i++) {
        double x = ((double)rand() / (double)RAND_MAX) * 20 - 10;
        Vector xi(1);
        xi(0) = (x - lower_bound_x1) / (upper_bound_x1 - lower_bound_x1);
        x_vec.push_back(xi);
        double yi;
        if (x < -2 && x >= -10) {
            yi = -2.186 * x - 12.864;
        } else if (x >= -2 && x < 0) {
            yi = 4.246 * x;
        } else if (x >= 0 && x <= 10) {
            yi = 10 * exp(-0.05 * x - 0.5) * sin((0.03*x + 0.7)*x);
        } else {
            std::cerr << "Error when generating data" << std::endl;
        }
        y.push_back((yi - lower_bound_y) / (upper_bound_y - lower_bound_y));
    }
}
