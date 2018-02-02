#include "WaveletFunction.h"
#include <stdexcept>
#include <math.h>

WaveletFunction::WaveletFunction() {

}

WaveletFunction::~WaveletFunction() {

}

double WaveletFunction::eval(const Vector& x) const {
    double sigma = 0.1;
    double res = 1;
    for (int i = 0; i<x.getN(); i++) {
        res *= -x(i) * std::exp(-std::pow(x(i), 2)/2.0);
        // res *= 2.0 / (std::sqrt(3.0 * sigma) * std::pow(M_PI, 1.0 / 4.0)) * (1 - x(i) * x(i) / (sigma * sigma)) * std::exp(- (x(i) * x(i)) / (2.0 * sigma * sigma)); 
    }
    return res;
}

Vector WaveletFunction::evalDev(const Vector& x) const {
    double sigma = 0.1;
    Vector gradient(x.getN());
    for (int i = 0; i<x.getN(); i++) {
        if (x(i) != 0) {
            gradient(i) = -(x(i) - 1.0/x(i)) * eval(x);
        } else {
            throw std::logic_error("Invalid wavelet derivative");
        }
        // gradient(i) = 2.0 / (std::sqrt(3.0 * sigma) * std::pow(M_PI, 1.0 / 4.0)) * (- 3.0 * x(i) / (sigma * sigma) + pow(x(i), 3) / pow(sigma, 4)) * std::exp(- (x(i) * x(i)) / (2.0 * sigma * sigma)); 
    }
    return gradient;
}
