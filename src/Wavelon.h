#ifndef WAVELON
#define WAVELON

#include <vector>
#include "Vector.h"
#include "Matrix.h"
#include "WaveletFunction.h"

class Wavelon
{
public:
    Wavelon(int input_dim);
    ~Wavelon();

    double computeOutput(const Vector& x);
    void updateParameters(const Vector& x, const double error);

protected:

    Matrix m_D;
    Matrix m_R;
    Vector m_t;
    double m_w;
    WaveletFunction m_phi;

    double m_z;

};

#endif
