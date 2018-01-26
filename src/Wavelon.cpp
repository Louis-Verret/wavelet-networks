#include "WaveletNetwork.h"

Wavelon::Wavelon(int input_dim) :
 m_D(input_dim, input_dim),
 m_t(input_dim),
 m_w(0),
 m_phi()
{
    //m_t.fillRandomly();
    //m_D.fillRandomly();
    m_R = Matrix::identity(input_dim, input_dim);
}

Wavelon::~Wavelon() {

}

double Wavelon::computeOutput(const Vector& x) {
    m_z = m_D * m_R * (x - m_t);
    return m_w * m_phi.eval(m_z);
}

void Wavelon::updateParameters(const Vector& x, const double error, const double learning_rate) {
    double delta_w = learning_rate * error * m_phi.eval(m_z);
    Vector delta_t = -error * m_w * m_R.transpose() * m_D * m_phi.evalDev(m_z);
    Vector delta_s = -error * m_w * m_D * m_D * Matrix::diag(m_R * (x - m_t)) * m_phi.evalDev(m_z);
    Matrix delta_R = error * m_w * m_D * m_phi.evalDev(m_z).mult(x - m_t);

    m_w = m_w - learning_rate * delta_w;
    m_t = m_t - learning_rate * delta_t;
    //std::cout << m_D << std::endl;
    m_D = (m_D.invDiag() - learning_rate * Matrix::diag(delta_s)).invDiag();
    m_R = m_R - learning_rate * delta_R;
    // std::cout << m_R << std::endl;
    //m_R = m_R.orthogonalization();
    // std::cout << "After projection " << m_R << std::endl;
}
