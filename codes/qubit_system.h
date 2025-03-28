#ifndef QUBIT_SYSTEM_H
#define QUBIT_SYSTEM_H

#include <complex>
#include"matrice.h"

class qubit;


class qubit_system {
public:
    virtual void display() const = 0;
    virtual std::complex<double> operator|(const qubit_system& q2) const = 0;  // Produit scalaire
    virtual matrice operator&(const qubit_system& q2) const = 0;    // Produit tensoriel
};

#endif
