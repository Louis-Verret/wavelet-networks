#include "WaveletNetwork.h"
#include "Utils.h"



WaveletNetwork::WaveletNetwork() : m_bar_g(0) {
}

WaveletNetwork::~WaveletNetwork() {
    for (int i = 0; i<m_nb_wavelons; i++) {
        delete m_wavelons[i];
    }
}

WaveletNetwork::WaveletNetwork(int input_dim, int nb_wavelons) : m_bar_g(0) {
    m_nb_wavelons = nb_wavelons;
    for (int i = 0; i<nb_wavelons; i++) {
        Wavelon* wavelon = new Wavelon(input_dim);
        m_wavelons.push_back(wavelon);
    }
}

void WaveletNetwork::init(std::vector<Vector>& x, std::vector<double>& y) {
    int s = x.size();
    int n = x[0].getN();
    std::vector<double> xi(s);
    std::vector<double> sorted_y = y;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < s; j++) {
            xi[j] = x[j](i);
        }
        // insertion sort
        for (int j = 0; j < s; j++) {
            double x = xi[j];
            double y = sorted_y[j];
            int k = j;
            while (k > 0 and xi[k - 1] > x) {
                xi[k] = xi[k - 1];
                sorted_y[k] = sorted_y[k - 1];
                k = k - 1;
            }
            xi[k] = x;
            sorted_y[k] = y;
        }
        for (int j = 0; j < s; j++) {
        }
    }
    // auto p = sort_permutation(x.getCoefficients(), [](T const& a, T const& b){ return a < b });
}


void WaveletNetwork::fit(std::vector<Vector>& x, std::vector<double>& y, const double learning_rate, int epoch) {
    for (int t = 0; t<epoch; t++) {
        double mse = 0;
        for (unsigned int i = 0; i<x.size(); i++) {
            double g = propagate(x[i]);
            double error = g - y[i];
            backpropagate(x[i], error, learning_rate);
            mse += std::pow(error, 2)/2.0;
        }
        std::cout << "Epoch: " << t+1 << "/" << epoch << std::endl;
        std::cout << " Mean Error: " << mse/x.size() << std::endl;
    }
    //std::cout << "Predicted/Label: " << m_a.back() << " " << batches_y[nb_batches-1] << std::endl;
}

double WaveletNetwork::propagate(const Vector& x) {
    double g = 0;
    for (int i = 0; i<m_nb_wavelons; i++) {
        g += m_wavelons[i]->computeOutput(x);
    }
    return g + m_bar_g;
}

void WaveletNetwork::backpropagate(const Vector& x, const double error, const double learning_rate) {
    for (int i = 0; i<m_nb_wavelons; i++) {
        m_wavelons[i]->updateParameters(x, error, learning_rate);
    }
    m_bar_g -= learning_rate * error;
}
