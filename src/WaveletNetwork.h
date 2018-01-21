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

    void fit(std::vector<Vector>& x, std::vector<double>& y, int epoch);

protected:

    std::vector<Wavelon*> m_wavelons;
    double m_bar_g;
    int m_nb_wavelons;

    double propagate(const Vector& x);
    void backpropagate(const Vector& x, const double error);

};

#endif