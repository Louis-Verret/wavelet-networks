#include "WaveletFunction.h"
#include <stdexcept>

WaveletFunction::WaveletFunction() {

}

WaveletFunction::~WaveletFunction() {

}

double WaveletFunction::eval(const Vector& x) const {
    double res = 1;
    for (int i = 0; i<x.getN(); i++) {
        res *= -x(i) * std::exp(-std::pow(x(i), 2)/2.0);
    }
    return res;
}

Vector WaveletFunction::evalDev(const Vector& x) const {
    Vector gradient(x.getN());
    for (int i = 0; i<x.getN(); i++) {
        if (x(i) != 0) {
            gradient(i) = -(x(i) - 1.0/x(i)) * eval(x);
        } else {
            throw std::logic_error("Invalid wavelet derivative");
        }
    }
    return gradient;
}
