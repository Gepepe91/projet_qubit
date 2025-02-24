#ifndef QUBIT_H
#define QUBIT_H

#include <complex>
#include <cmath>

using namespace std;

class Qubit {
private:
    complex<double> psi0;  // Amplitude pour |0>
    complex<double> psi1;  // Amplitude pour |1>
    double theta;          // Paramètre theta
    double phi;            // Paramètre phi

public:
    // Constructeur pour initialiser un qubit avec des paramètres θ et φ
    Qubit(double theta, double phi);

    // Getter pour obtenir l'amplitude de |0>
    complex<double> getPsi0() const;

    // Getter pour obtenir l'amplitude de |1>
    complex<double> getPsi1() const;

    // Setter pour modifier θ et φ
    void setThetaPhi(double theta, double phi);

    // Calcul des probabilités pour |0> et |1>
    double prob_psi0() const;
    double prob_psi1() const;
};

#endif