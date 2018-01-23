#include "WaveletNetwork.h"
#include "Utils.h"



WaveletNetwork::WaveletNetwork() : m_bar_g(0) {
}

WaveletNetwork::~WaveletNetwork() {
    for (int i = 0; i<m_nb_wavelons; i++) {
        delete m_wavelons[i];
    }
}

WaveletNetwork::WaveletNetwork(int input_dim, int nb_wavelons) : m_bar_g(0), m_nb_wavelons(nb_wavelons) {
    for (int i = 0; i<nb_wavelons; i++) {
        Wavelon* wavelon = new Wavelon(input_dim);
        m_wavelons.push_back(wavelon);
    }
}

void WaveletNetwork::recursively_init(std::vector<double>& xi, std::vector<double>& y, Vector t_vec, Vector s_vec, int wavelon_index) {
    // insertion sort
    int s = xi.size();
    std::vector<double> sorted_y = y;
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
    // derivative computation
    std::vector<double> rhovar(s - 1);
    for (int j = 0; j < s - 1; j++) {
        rhovar[j] = std::abs((sorted_y[j + 1] - sorted_y[j]) / (xi[j + 1] - xi[j]));
        // std::cout << xi[j] << std::endl;
        // std::cout << rhovar[j] << std::endl;
    }
    //integral computation
    double integral = 0;
    for (int j = 0; j < s - 1; j++) {
        integral += (xi[j + 1] - xi[j]) * rhovar[j];  
        // std::cout << integral << std::endl;
    }
    //rhovar computation
    std::vector<double> rho(s - 1);
    for (int j = 0; j < s - 1; j++) {
        rho[j] = rhovar[j] / integral;
        // std::cout << rho[j] << std::endl;
    }
    //computing p
    double p = 0;
    for (int j = 0; j < s - 1; j++) {
        p += (xi[j + 1] - xi[j])  * (xi[j] * rho[j]);
        //std::cout << p << std::endl;
    }
    // parameters initialization
    //t_vec(wavelon_index) = p;
    //s_vec(wavelon_index) = 0.5 * (xi[s - 1] - xi[0]);

    // spliting data in two
    std::vector<double> left_xi; std::vector<double> right_xi;
    std::vector<double> left_y; std::vector<double> right_y;
    for (int j = 0; j < s; j++) {
        if (xi[j] < p) {
            left_xi.push_back(xi[j]);
            left_y.push_back(y[j]);
        } else {
            right_xi.push_back(xi[j]);
            right_y.push_back(y[j]);
        }
    }
    // std::cout << left_xi.size() << std::endl;
    // std::cout << right_xi.size() << std::endl;
    
    if (wavelon_index < m_nb_wavelons) {
        recursively_init(left_xi, left_y, t_vec, s_vec, wavelon_index + 1);
        recursively_init(right_xi, right_y, t_vec, s_vec, wavelon_index + 1);
    }
}

void WaveletNetwork::init(std::vector<Vector>& x, std::vector<double>& y) {
    int s = x.size();
    int n = x[0].getN();
    std::vector<double> xi(s);
    for (int i = 0; i < n; i++) {
        // get vector of i-nth coordinates of x
        for (int j = 0; j < s; j++) {
            xi[j] = x[j](i);
        }
        Vector t_vec(m_nb_wavelons);
        Vector s_vec(m_nb_wavelons);
        recursively_init(xi, y, t_vec, s_vec, 0);
    }
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
