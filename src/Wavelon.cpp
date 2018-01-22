#include "WaveletNetwork.h"

Wavelon::Wavelon(int input_dim) :
 m_D(input_dim, input_dim),
 m_R(input_dim, input_dim),
 m_t(input_dim),
 m_w(0.01),
 m_phi()
{
    m_D.fillRandomly();
    m_R.fillRandomly();
    m_t.fillRandomly();
}

Wavelon::~Wavelon() {

}

double Wavelon::computeOutput(const Vector& x) {
    m_z = m_w * (m_D * m_R * (x - m_t));
    return m_phi.eval(m_z);
}

void Wavelon::updateParameters(const Vector& x, const double error, const double learning_rate) {
    double delta_w = learning_rate * error * m_phi.eval(m_z);
    Vector delta_t = -error * m_w * m_R.transpose() * m_D * m_phi.evalDev(m_z);
    Vector delta_s = -error * m_w * m_D * m_D * Matrix::diag(m_R * (x - m_t)) * m_phi.evalDev(m_z);
    Matrix delta_R = error * m_w * m_D * m_phi.evalDev(m_z).mult(x - m_t);

    m_w = m_w - learning_rate * delta_w;
    m_t = m_t - learning_rate * delta_t;
    m_D = m_D - learning_rate * Matrix::diag(delta_s.inv());
    m_R = m_R - learning_rate * delta_R;
}
