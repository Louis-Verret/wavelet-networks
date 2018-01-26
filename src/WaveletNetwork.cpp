#include "WaveletNetwork.h"
#include "Utils.h"
#include <queue>

WaveletNetwork::WaveletNetwork() : m_bar_g(0) {
}

WaveletNetwork::~WaveletNetwork() {
    for (int i = 0; i<m_nb_wavelons; i++) {
        delete m_wavelons[i];
    }
}

WaveletNetwork::WaveletNetwork(int input_dim, int nb_wavelons) : m_bar_g(0), m_nb_wavelons(nb_wavelons) {
    for (int i = 0; i<nb_wavelons; i++) {
        Wavelon* wavelon = new Wavelon(input_dim);
        m_wavelons.push_back(wavelon);
    }
}

void WaveletNetwork::recursively_init(std::vector<double>& xi, std::vector<double>& y, std::vector<double>& vec_s, std::vector<double>& vec_t) {
    std::queue<std::vector<double> > q_xi;
    std::queue<std::vector<double> > q_y;
    q_xi.push(xi);
    q_y.push(y);
    int l = 0;
    while (q_xi.size() > 0) {
        std::vector<double> sub_xi = q_xi.front();
        std::vector<double> sub_y = q_y.front();
        q_xi.pop(); q_y.pop();
        int s = sub_xi.size();
        //std::cout << l << " " << s << std::endl;
        // derivative computation
        std::vector<double> rhovar(s - 1);
        for (int j = 0; j < s - 1; j++) {
            rhovar[j] = std::abs((sub_y[j + 1] - sub_y[j]) / (sub_xi[j + 1] - sub_xi[j]));
            // std::cout << xi[j] << std::endl;
            // std::cout << sub_xi[j] << std::endl;
            // std::cout << rhovar[j] << std::endl;
        }
        //integral computation
        double integral = 0;
        for (int j = 0; j < s - 1; j++) {
            integral += (sub_xi[j + 1] - sub_xi[j]) * rhovar[j];
            // std::cout << integral << std::endl;
        }
        //rhovar computation
        std::vector<double> rho(s - 1);
        for (int j = 0; j < s - 1; j++) {
            rho[j] = rhovar[j] / integral;
            // std::cout << rho[j] << std::endl;
        }
        //computing p
        double p = 0;
        for (int j = 0; j < s - 1; j++) {
            p += (sub_xi[j + 1] - sub_xi[j])  * (sub_xi[j] * rho[j]);
            // std::cout << p << std::endl;
        }
        if (l < m_nb_wavelons) {
            vec_t[l] = p;
            vec_s[l] = 0.5 * (sub_xi[s - 1] - sub_xi[0]);
            // spliting data in two
            std::vector<double> left_xi; std::vector<double> right_xi;
            std::vector<double> left_y; std::vector<double> right_y;
            for (int j = 0; j < s; j++) {
                if (sub_xi[j] < p) {
                    left_xi.push_back(sub_xi[j]);
                    left_y.push_back(sub_y[j]);
                } else {
                    right_xi.push_back(sub_xi[j]);
                    right_y.push_back(sub_y[j]);
                }
            }
            q_xi.push(left_xi); q_xi.push(right_xi);
            q_y.push(left_y); q_y.push(right_y);
            l++;
        } else {
            break;
        }

    }

}

void WaveletNetwork::init(std::vector<Vector>& x, std::vector<double>& y) {
    int s = x.size();
    int n = x[0].getN();
    std::vector<double> xi(s);
    std::vector<std::vector<double> > all_t;
    std::vector<std::vector<double> > all_s;
    for (int i = 0; i < n; i++) {
        // get vector of i-nth coordinates of x
        for (int j = 0; j < s; j++) {
            xi[j] = x[j](i);
        }
        // insertion sort
        std::vector<double> sorted_y = y;
        for (int j = 0; j < s; j++) {
            double x = xi[j];
            double y = sorted_y[j];
            int k = j;
            while (k > 0 and xi[k - 1] > x) {
                xi[k] = xi[k - 1];
                sorted_y[k] = sorted_y[k - 1];
                k = k - 1;
            }
            xi[k] = x;
            sorted_y[k] = y;
        }
        std::vector<double> t_vec(m_nb_wavelons);
        std::vector<double> s_vec(m_nb_wavelons);
        recursively_init(xi, y, t_vec, s_vec);
        all_t.push_back(t_vec);
        all_s.push_back(s_vec);
    }
    for (int l = 0; l < m_nb_wavelons; l++) {
        Vector t(n);
        Matrix D(n, n);
        D.fillWithZero();
        // std::cout << "Wavelon number " << l << std::endl;
        for (int i = 0; i < n; i++) {
            t(i) = all_t[i][l];
            // std::cout << "Coordinate "  << l << " --> " << t(i) << std::endl;
            D(i, i) = 1 / all_s[i][l];
            // std::cout << D(i, i) << std::endl;
        }
        m_wavelons[l]->setT(t);
        m_wavelons[l]->setD(D);
    }
    // gbar computation
    double g_bar = 0;
    for (int j = 0; j < s; j++) {
        g_bar += y[j];
    }
    // std::cout << g_bar / s << std::endl;
    m_bar_g = g_bar / s;
}


void WaveletNetwork::fit(std::vector<Vector>& x, std::vector<double>& y, const double learning_rate, int epoch) {
    double g;
    for (int t = 0; t<epoch; t++) {
        double mse = 0;
        for (unsigned int i = 0; i<x.size(); i++) {
            g = propagate(x[i]);
            double error = g - y[i];
            backpropagate(x[i], error, learning_rate);
            mse += std::pow(error, 2)/2.0;
        }
        std::cout << "Epoch: " << t+1 << "/" << epoch << std::endl;
        std::cout << " Mean Error: " << mse/x.size() << std::endl;
    }
    std::cout << "Predicted/Label: " << g << " " << y[y.size()-1] << std::endl;
}

double WaveletNetwork::propagate(const Vector& x) {
    double g = 0;
    for (int i = 0; i<m_nb_wavelons; i++) {
        g += m_wavelons[i]->computeOutput(x);
    }
    return g + m_bar_g;
}

void WaveletNetwork::backpropagate(const Vector& x, const double error, const double learning_rate) {
    for (int i = 0; i<m_nb_wavelons; i++) {
        m_wavelons[i]->updateParameters(x, error, learning_rate);
    }
    m_bar_g -= learning_rate * error;
}
