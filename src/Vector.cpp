#include "Vector.h"
#include <iostream>
#include <cmath>
#include <stdexcept>

Vector::Vector() : m_n(0) {
    //m_coefficients = new double[0];
    m_coefficients = std::vector<double>(0);
}

Vector::~Vector() {
    //delete[] m_coefficients;
}

Vector::Vector(int n) :
        m_n(n)
{
    m_coefficients = std::vector<double>(n);
}

Vector::Vector(int n, double val) :
        m_n(n)
{
    m_coefficients = std::vector<double>(n, val);
    // m_coefficients = new double[n];
    // #pragma omp parallel shared(val, n)
    // {
    //     #pragma omp for
    //     for (int i = 0; i < n; i++) {
    //         m_coefficients[i] = val;
    //     }
    // }
}

// Vector::Vector(const Vector& vec) {
//     int n = vec.getN();
//     m_coefficients = new double[n];
//     std::copy(vec.m_coefficients, vec.m_coefficients + n, m_coefficients);
// }

const double &Vector::operator()(int i) const {
    if (i < m_n)
        return m_coefficients[i];
    throw std::logic_error("Invalid vector element");
}

double &Vector::operator()(int i) {
    if (i < m_n)
        return m_coefficients[i];
    throw std::logic_error("Invalid vector element");
}

// Vector& Vector::operator=(const Vector& vec) {
//     int n = vec.getN();
//     if (n != m_n) {
//         delete[] m_coefficients;
//         m_coefficients = new double[n];
//     }
//     std::copy(vec.m_coefficients, vec.m_coefficients + n, m_coefficients);
//     return *this;
// }

Vector Vector::operator+(const Vector &v) const {
    if (v.getN() != m_n) {
        throw std::logic_error("Invalid size for vector addition");
    }
    Vector res(m_n);
    #pragma omp parallel shared(res, v)
    {
        #pragma omp for
        for (int i = 0; i < m_n; i++) {
            res(i) = m_coefficients[i] + v(i);
        }
    };
    return res;
}

Vector Vector::operator-(const Vector &v) const {
    if (v.getN() != m_n) {
        throw std::logic_error("Invalid size for vector substraction");
    }
    Vector res(m_n);
    #pragma omp parallel shared(res, v)
    {
        #pragma omp for
        for (int i = 0; i < m_n; i++) {
            res(i) = m_coefficients[i] - v(i);
        }
    };
    return res;
}

Vector Vector::operator*(const Vector &v) const {
    if (v.getN() != m_n) {
        throw std::logic_error("Invalid size for Hadamard vector product");
    }
    Vector res(m_n);
    #pragma omp parallel shared(res, v)
    {
        #pragma omp for
        for (int i = 0; i < m_n; i++) {
            res(i) = m_coefficients[i] * v(i);
        }
    };
    return res;
}

Vector Vector::operator/(const Vector &v) const {
    if (v.getN() != m_n) {
        throw std::logic_error("Invalid size for vector division");
    }
    Vector res(m_n);
    #pragma omp parallel shared(res, v)
    {
        #pragma omp for
        for (int i = 0; i < v.getN(); i++) {
            res(i) = m_coefficients[i] / v(i);
        }
    };
    return res;
}

Vector Vector::operator/(const double coeff) const {
    Vector res(m_n);
    #pragma omp parallel shared(res)
    {
        #pragma omp for
        for (int i = 0; i < m_n; i++) {
            res(i) = m_coefficients[i] / coeff;
        }
    };
    return res;
}

Vector Vector::operator+(const double coeff) const {
    Vector res(m_n);
    #pragma omp parallel shared(res)
    {
        #pragma omp for
        for (int i = 0; i < m_n; i++) {
            res(i) = m_coefficients[i] + coeff;
        }
    };
    return res;
}

void Vector::fillRandomly() {
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < m_n; i++) {
            double r = ((double) rand() / (double) RAND_MAX);
            m_coefficients[i] = r;
        }
    };
}

void Vector::fillWithZero() {
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < m_n; i++) {
            m_coefficients[i] = 0;
        }
    }
}

Vector Vector::sqrt() const {
    Vector res(m_n);
    #pragma omp parallel shared(res)
    {
        #pragma omp for
        for (int i = 0; i < m_n; i++) {
            if (m_coefficients[i] < 0) {
                throw std::logic_error("Invalid Vector sqrt");
            }
            res(i) = std::sqrt(m_coefficients[i]);
        }
    };
    return res;
}

Vector Vector::inv() const {
    Vector res(m_n);
    #pragma omp parallel shared(res)
    {
        #pragma omp for
        for (int i = 0; i < m_n; i++) {
            if (m_coefficients[i] == 0) {
                throw std::logic_error("Invalid Vector inv");
            }
            res(i) = 1.0/m_coefficients[i];
        }
    };
    return res;
}

Matrix Vector::mult(const Vector& vec) const {
    Matrix result(m_n, vec.getN());
    int i, j;
    #pragma omp parallel shared(result) private(i, j)
    {
        #pragma omp for
        for (i = 0; i < m_n; i++) {
            for (j = 0; j < m_n; j++) {
                result(i, j) = m_coefficients[i] * vec(j);
            }
        }
    }
    return result;
}

Vector operator*(const double coeff, const Vector& v) {
    Vector res(v.getN());
    #pragma omp parallel shared(res, v)
    {
        #pragma omp for
        for (int i = 0; i < v.getN(); i++) {
            res(i) = coeff * v(i);
        }
    };
    return res;
}

std::ostream& operator << (std::ostream& out, const Vector& v) {
    int n = v.getN();
    for (int i = 0; i < n; i++) {
        out << v(i) << " ";
    }
    return out;
}
