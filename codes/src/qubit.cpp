#include "qubit.h"
#include <complex>
#include <iostream>
#include"matrice.h"

using complexe = std::complex<double>;

// Constructors
qubit::qubit() : theta(0), phi(0) {}

qubit::qubit(complexe alpha_, complexe beta_) : alpha(alpha_), beta(beta_) {
    synchronize_alpha_beta_to_theta_phi();
}

qubit::qubit(double theta_, double phi_) : theta(theta_), phi(phi_) {
    synchronize_theta_phi_to_alpha_beta();
}

// Synchronization between alpha/beta and theta/phi
void qubit::synchronize_alpha_beta_to_theta_phi() {
    if (std::abs(alpha) >= 1 - 1e-12) {
        theta = 0.0;
        phi = 0.0;
    } else if (std::abs(beta) >= 1 - 1e-12) {
        theta = M_PI;
        phi = 0.0;
    } else {
        theta = 2 * std::acos(std::abs(alpha));
    }
    phi = std::arg(beta) - std::arg(alpha);
}

void qubit::synchronize_theta_phi_to_alpha_beta() {
    alpha = std::cos(theta / 2);
    beta = std::exp(complexe(0, phi)) * std::sin(theta / 2);
}

// Getters
complexe qubit::get_alpha() const { return alpha; }
complexe qubit::get_beta() const { return beta; }
double qubit::get_theta() const { return theta; }
double qubit::get_phi() const { return phi; }
double qubit::get_abs_alpha2() const { return std::norm(alpha); }
double qubit::get_abs_beta2() const { return std::norm(beta); }

// Setters
void qubit::set_alpha(complexe alpha_) { alpha = alpha_; synchronize_alpha_beta_to_theta_phi(); }
void qubit::set_beta(complexe beta_) { beta = beta_; synchronize_alpha_beta_to_theta_phi(); }
void qubit::set_theta(double theta_) { theta = theta_; synchronize_theta_phi_to_alpha_beta(); }
void qubit::set_phi(double phi_) { phi = phi_; synchronize_theta_phi_to_alpha_beta(); }

// Normalization
void qubit::normalize() {
    double norm_const = std::sqrt(std::norm(alpha) + std::norm(beta));
    alpha /= norm_const;
    beta /= norm_const;
}

// Display method
void qubit::display() const {
    std::cout << std::endl;
    std::cout << alpha << " coefficient of |0>" << std::endl;
    std::cout << beta << " coefficient of |1>" << std::endl;
    std::cout << std::endl;
}

void qubit::display_angles() {
    std::cout << "Theta: " << theta << " | Phi: " << phi << std::endl;
}

// Operator overloads
qubit qubit::operator*(const complexe& scalaire) const {
    return qubit(alpha * scalaire, beta * scalaire);
}

qubit qubit::operator*(const double& scalaire) const{
    return qubit(alpha * scalaire, beta * scalaire);
}

qubit qubit::operator+(const qubit& r) const{
    return qubit(alpha + r.alpha, beta + r.beta);
}

std::complex<double> qubit::operator|(const qubit_system& q2) const {
    const qubit* q = dynamic_cast<const qubit*>(&q2);
    if (!q) throw std::runtime_error("Invalid qubit_system conversion to qubit.");
    return std::conj(alpha) * q->alpha + std::conj(beta) * q->beta;
}

matrice qubit::operator&(const qubit_system& q2) const {
    const qubit* q = dynamic_cast<const qubit*>(&q2);
    if (!q) throw std::runtime_error("Invalid qubit_system conversion to qubit.");
    return matrice(alpha * std::conj(q->alpha), alpha * std::conj(q->beta),
                   beta * std::conj(q->alpha), beta * std::conj(q->beta));
}
