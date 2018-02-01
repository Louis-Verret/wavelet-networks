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
    std::vector<double> y_predict_1(y_train.size(), 0);

    int nb_wavelons = 5;
    int input_dim = 1;
    WaveletNetwork net(input_dim, nb_wavelons);

    std::cout << "Initializing the network" << std::endl;
    net.init(x_train, y_train);

    std::cout << "Fitting the data" << std::endl;
    net.fit(x_train, y_train, 0.05, 10000);

    std::cout << "Evaluating the data" << std::endl;
    net.evaluate(x_test, y_predict);

    std::ofstream file("../data/prediction_func1_test_1.data");
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

    std::cout << "Evaluating the data" << std::endl;
    net.evaluate(x_train, y_predict_1);

    std::ofstream file1("../data/prediction_func1_train_1.data");
    if (file1.is_open())
    {
        int s = x_train.size();
        int n = x_train[0].getN();
        for (int i = 0; i < s; i++) {
            for (int j = 0; j < n; j++) {
                file1 << x_train[i](j) << " ";
            }
            file1 << y_predict_1[i] << " ";
            file1 << y_train[i] << std::endl;
        }
        file1.close();
    }
    else std::cout << "Unable to open file";
    return 0;

    // Matrix A(2,2);
    // A.fillRandomly();
    // // for (int i = 0; i < 5; i++) {
    // //     for (int j = 0; j < 5; j++) {
    // //         A(i,j) = i - j;
    // //     }
    // // }
    // std::cout << A << std::endl;
    // std::cout << A.transpose() * A << std::endl;
    // Matrix B = A.orthogonalization();
    // std::cout << B.transpose() * B  << std::endl;
    // return 0;
}
