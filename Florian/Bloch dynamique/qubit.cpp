#include "qubit.h"
using namespace std;

// Constructeur qui initialise l'état du qubit avec θ et φ
Qubit::Qubit(double theta, double phi) : theta(theta), phi(phi) {
    // Initialisation de l'état |ψ⟩
    psi0 = cos(theta / 2);  // Amplitude de |0>
    psi1 = exp(complex<double>(0, phi)) * sin(theta / 2);  // Amplitude de |1>
}

// Getter pour obtenir l'amplitude de |0>
complex<double> Qubit::getPsi0() const {
    return psi0;
}

// Getter pour obtenir l'amplitude de |1>
complex<double> Qubit::getPsi1() const {
    return psi1;
}

// Setter pour modifier θ et φ
void Qubit::setThetaPhi(double theta, double phi) {
    this->theta = theta;
    this->phi = phi;
    // Recalculer les amplitudes après modification des paramètres
    psi0 = cos(theta / 2);
    psi1 = exp(complex<double>(0, phi)) * sin(theta / 2);
}

// Calcul des probabilités pour |0> et |1>
double Qubit::prob_psi0() const {
    return norm(psi0);
}

double Qubit::prob_psi1() const {
    return norm(psi1);
}
