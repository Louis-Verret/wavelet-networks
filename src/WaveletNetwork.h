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

    /* Initialise the network parameters using the algorithm proposed in the paper */
    void init(std::vector<Vector>& x, std::vector<double>& y);

    /* Implements the learning algorithm using input vector x,
        known labels y, at a given learning rate and for a given number of iterations */
    void fit(std::vector<Vector>& x, std::vector<double>& y, const double learning_rate, int epoch);

    /* Outputs on y the value of the predicted value for the input x */
    void evaluate(std::vector<Vector>& x, std::vector<double>& y);

protected:
    std::vector<Wavelon*> m_wavelons;     /* collection of wavelons of the model */
    double m_bar_g;                       /* additional parameters introduced in the paper */
    int m_nb_wavelons;                    /* number of wavelons /*

    /* Return the value computed after propagating the vector x through the network */
    double propagate(const Vector& x);

    /* Uses the given error to correct the parameters using stochastic gradient descent */
    void backpropagate(const Vector& x, const double error, const double learning_rate);

    /* Actual computation of the initialisation of the parameters of the network */
    void recursively_init(std::vector<double>& xi, std::vector<double>& y, std::vector<double>& vec_s, std::vector<double>& vec_t);

};

#endif
