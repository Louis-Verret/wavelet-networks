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
    void updateParameters(const Vector& x, const double error, const double learning_rate);

    Vector getT() const { return m_t; };
    void setT(const Vector& t) { m_t = t; };

protected:

    Matrix m_D;
    Matrix m_R;
    Vector m_t;
    double m_w;
    WaveletFunction m_phi;

    Vector m_z;

};

#endif
