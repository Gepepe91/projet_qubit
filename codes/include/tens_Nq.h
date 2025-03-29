#ifndef TENS_NQ_H
#define TENS_NQ_H

#include "qubit_system.h"
#include <vector>
#include <complex>

class tens_Nq : public qubit_system {
private:
    int n;  // Nombre de qubits
    std::vector<std::complex<double>> amplitudes; // 2^n coefficients complexes

public:
    // Constructeur
    explicit tens_Nq(int n); // Initialise un état à n qubits

    // Redéfinition des méthodes virtuelles
    void display() const override;
    std::complex<double> operator|(const qubit_system& qs) const override;
    matrice operator&(const qubit_system& qs) const override;

    // Méthode de normalisation
    void normalize();

    // Getters
    int get_n() const { return n; }
    const std::vector<std::complex<double>>& get_amplitudes() const { return amplitudes; }
};

#endif
