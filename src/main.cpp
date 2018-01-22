#include "WaveletNetwork.h"
#include "Utils.h"

#include <iostream>


int main(int argc, char **argv) {

    std::vector<Vector> x;
    std::vector<double> y;

    generateData(x, y, 100);

    int nb_wavelons = 1;
    int input_dim = 1;
    WaveletNetwork net(nb_wavelons, input_dim);

    net.init(x, y);
    net.fit(x, y, 0.001, 2000);

    return 0;
}
