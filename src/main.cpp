#include "WaveletNetwork.h"
#include "Utils.h"

#include <iostream>
#include "Matrix.h"


int main(int argc, char **argv) {

    std::vector<Vector> x_train;
    std::vector<double> y_train;
    generateData(x_train, y_train, 100);

    // std::vector<Vector> x_test;
    // std::vector<double> y_test;
    // generateData(x_test, y_test, 1000);

    int nb_wavelons = 8;
    int input_dim = 1;
    WaveletNetwork net(input_dim, nb_wavelons);

    std::cout << "Initializing the network" << std::endl;
    net.init(x_train, y_train);

    std::cout << "Fitting the data" << std::endl;
    net.fit(x_train, y_train, 0.1, 1500);

    // Matrix A = 2 * Matrix::identity(5, 5);
    // // for (int i = 0; i < 5; i++) {
    // //     for (int j = 0; j < 5; j++) {
    // //         A(i,j) = i - j;
    // //     }
    // // }
    // std::cout << A << std::endl;
    // std::cout << A.orthogonalization() << std::endl;
    // return 0;
}
