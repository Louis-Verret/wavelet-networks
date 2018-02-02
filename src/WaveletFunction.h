#ifndef WAVELETFUNCTION
#define WAVELETFUNCTION

#include <vector>
#include <cmath>
#include "Vector.h"

class WaveletFunction
{
public:
    WaveletFunction();
    ~WaveletFunction();

    /* Evaluate the wavelet function on the vector x */
    double eval(const Vector& x) const;

    /* Evaluate the gradient of the wavelet function on the vector x */
    Vector evalDev(const Vector& x) const;

};

#endif
