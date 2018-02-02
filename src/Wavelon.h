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

    /* Compute the output of the wavelon for the input x */
    double computeOutput(const Vector& x);

    /* Updates parameters of the wavelon according to the error, using gradient descent computation */
    void updateParameters(const Vector& x, const double error, const double learning_rate);

    Vector getT() const { return m_t; };        /* Return the translation vector */
    void setT(const Vector& t) { m_t = t; };    /* Set the translation vector */
    Matrix getD() const { return m_D; };        /* Return the dilatation matrix */
    void setD(const Matrix& D) { m_D = D; };    /* Set the dilatation matrix */

protected:
    Matrix m_D;                 /* dilatation matrix */
    Matrix m_R;                 /* rotation matrix */
    Vector m_t;                 /* translation vector */
    double m_w;                 /* weight */
    WaveletFunction m_phi;      /* mother wavelet function */

    Vector m_z;

};

#endif
