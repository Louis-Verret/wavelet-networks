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

    double eval(const Vector& x) const;
    double evalDev(const Vector& x) const;

};

#endif