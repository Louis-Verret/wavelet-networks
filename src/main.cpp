#include "WaveletNetwork.h"
#include "Utils.h"

#include <iostream>
#include "Matrix.h"
#include <fstream>


int main(int argc, char **argv) {

    std::vector<Vector> x_train;
    std::vector<double> y_train;
    generateData3(x_train, y_train, 300);

    std::vector<Vector> x_test;
    std::vector<double> y_test;
    generateData3(x_test, y_test, 200);

    std::vector<double> y_predict(y_test.size(), 0);

    int nb_wavelons = 11;
    int input_dim = 1;
    double learning_rate = 0.05;
    WaveletNetwork net(input_dim, nb_wavelons);

    std::cout << "Initializing the network" << std::endl;
    net.init(x_train, y_train);

    std::cout << "Fitting the data" << std::endl;
    net.fit(x_train, y_train, learning_rate, 10000);

    std::cout << "Evaluating the data" << std::endl;
    net.evaluate(x_test, y_predict);

    std::ofstream file("../data/prediction_func1_test.data");
    if (file.is_open())
    {
        int s = x_test.size();
        int n = x_test[0].getN();
        for (int i = 0; i < s; i++) {
            for (int j = 0; j < n; j++) {
                file << x_test[i](j) << " ";
            }
            file << y_predict[i] << " ";
            file << y_test[i] << std::endl;
        }
        file.close();
    }
    else std::cout << "Unable to open file";
}
