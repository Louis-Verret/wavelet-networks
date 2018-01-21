#include "WaveletFunction.h"

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

double WaveletFunction::evalDev(const Vector& x) const {

}
