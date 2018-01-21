#include "WaveletNetwork.h"

Wavelon::Wavelon(int input_dim) :
 m_D(input_dim, input_dim),
 m_R(input_dim, input_dim),
 m_t(input_dim),
 m_w(0),
 m_phi()
{
    m_D.fillRandomly();
    m_R.fillRandomly();
    m_t.fillRandomly();
}

Wavelon::~Wavelon() {

}

double Wavelon::computeOutput(const Vector& x) {
    m_z = m_phi.eval(m_w * (m_D * m_R * (x - m_t)));
    return m_z;
}

void Wavelon::updateParameters(const Vector& x, const double error) {

}
