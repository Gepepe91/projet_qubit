#ifndef qubit_H
#define qubit_H

#include "qubit_system.h"
#include <complex>

using complexe = std::complex<double>;

class qubit : public qubit_system {
private:
    complexe alpha, beta;
    double theta, phi;

    void synchr_alpha_beta_to_theta_phi();
    void synchr_theta_phi_to_alpha_beta();

    void normalize();

public:
    // Constructeurs
    qubit();
    qubit(std::complex<double> alpha_, std::complex<double> beta_);
    qubit(double theta_, double phi_);
    
    // Redéfinition des méthodes virtuelles
    void display() const override;
    std::complex<double> operator|(const qubit_system& q2) const override;
    matrice operator&(const qubit_system& q2) const override;
    void display_angles();

    //getters
    complexe get_alpha() const;
    complexe get_beta() const;
    double get_theta() const;
    double get_phi() const;
    double get_abs_alpha2() const;
    double get_abs_beta2() const;

    void set_alpha(complexe alpha_);
    void set_beta(complexe beta_);
    void set_theta(double theta_);
    void set_phi(double phi_);

    qubit operator*(const complexe& scalaire);
    qubit operator*(const double& scalaire);
    qubit operator+(const qubit& r);
    complexe operator|(const qubit& q2);
};

#endif
