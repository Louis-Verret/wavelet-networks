#include "WaveletNetwork.h"
#include "Utils.h"

#include <iostream>


int main(int argc, char **argv) {

    std::vector<Vector> x;
    std::vector<double> y;

    generateData2(x, y, 1000);

    int nb_wavelons = 4;
    int input_dim = 1;
    WaveletNetwork net(input_dim, nb_wavelons);

    net.init(x, y);
    net.fit(x, y, 5);

    return 0;
}
