#include "WaveletNetwork.h"
#include "Utils.h"

#include <iostream>


int main(int argc, char **argv) {

    std::vector<Vector> x;
    std::vector<double> y;

    generateData(x, y, 100);

    int nb_wavelons = 8;
    int input_dim = 1;
    WaveletNetwork net(input_dim, nb_wavelons);

    net.init(x, y);
    net.fit(x, y, 0.1, 1500);

    return 0;
}
