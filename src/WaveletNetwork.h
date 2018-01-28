#ifndef WAVELET_NETWORK
#define WAVELET_NETWORK

#include <vector>
#include "Vector.h"
#include "Matrix.h"
#include "Wavelon.h"
#include <cmath>
#include <iostream>

class WaveletNetwork
{
public:
    WaveletNetwork();
    WaveletNetwork(int nb_wavelons, int input_dim);
    ~WaveletNetwork();

    void init(std::vector<Vector>& x, std::vector<double>& y);
    void fit(std::vector<Vector>& x, std::vector<double>& y, const double learning_rate, int epoch);
    void evaluate(std::vector<Vector>& x, std::vector<double>& y);

protected:

    std::vector<Wavelon*> m_wavelons;
    double m_bar_g;
    int m_nb_wavelons;

    double propagate(const Vector& x);
    void backpropagate(const Vector& x, const double error, const double learning_rate);

    void recursively_init(std::vector<double>& xi, std::vector<double>& y, std::vector<double>& vec_s, std::vector<double>& vec_t);

};

#endif
