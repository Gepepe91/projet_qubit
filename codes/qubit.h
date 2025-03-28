#ifndef QUBIT_H
#define QUBIT_H

#include <complex>
#include"qubit_system.h"

using complexe = std::complex<double>;

class qubit: public qubit_system {
private:
    complexe alpha, beta;
    double theta, phi;

    void synchronize_alpha_beta_to_theta_phi();
    void synchronize_theta_phi_to_alpha_beta();

public:
    void normalize();

    // Constructors
    qubit();
    qubit(complexe alpha_, complexe beta_);
    qubit(double theta_, double phi_);

    // Virtual method overrides
    void display() const override;
    std::complex<double> operator|(const qubit_system& q2) const override;
    matrice operator&(const qubit_system& q2) const override;
    void display_angles();

    // Getters
    complexe get_alpha() const;
    complexe get_beta() const;
    double get_theta() const;
    double get_phi() const;
    double get_abs_alpha2() const;
    double get_abs_beta2() const;

    // Setters
    void set_alpha(complexe alpha_);
    void set_beta(complexe beta_);
    void set_theta(double theta_);
    void set_phi(double phi_);

    // Operator overloads
    qubit operator*(const complexe& scalaire) const;
    qubit operator*(const double& scalaire) const;
    qubit operator+(const qubit& r) const;
};

#endif
