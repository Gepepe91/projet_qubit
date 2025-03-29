#ifndef TENS_2Q_H
#define TENS_2Q_H

#include "qubit_system.h"
#include "qubit.h"

class tens_2q : public qubit_system {
private:
    qubit q1, q2; // Les deux qubits du système

public:
    // Constructeurs
    tens_2q();
    tens_2q(const qubit& q1_, const qubit& q2_);
    tens_2q(std::complex<double> alpha00, std::complex<double> alpha01, 
            std::complex<double> alpha10, std::complex<double> alpha11);

    // Redéfinition des méthodes virtuelles
    void display() const override;  
    complexe operator|(const qubit_system& q) const override;  // Produit scalaire
    matrice operator&(const qubit_system& q) const override;   // Produit tensoriel

    // Accesseurs
    qubit get_q1() const ;
    qubit get_q2() const ;
};

#endif
