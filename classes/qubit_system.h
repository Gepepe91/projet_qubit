#ifndef QUBIT_SYSTEM_H
#define QUBIT_SYSTEM_H

#include "matrice.h"
#include <complex>
#include <iostream>

using complexe = std::complex<double>;

class qubit_system {
public:
    virtual ~qubit_system() = default;  // Classe de base avec destructeur virtuel

    // Méthodes virtuelles pures à redéfinir dans les classes dérivées
    virtual void display() const = 0;
    virtual complexe operator|(const qubit_system& q2) const = 0;  // Produit scalaire
    virtual matrice operator&(const qubit_system& q2) const = 0;    // Produit tensoriel
};

#endif
