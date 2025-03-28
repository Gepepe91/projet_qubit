#include "tens_2q.h"
#include <iostream>

// Constructeurs
tens_2q::tens_2q() : q1(), q2() {}

tens_2q::tens_2q(const qubit& q1_, const qubit& q2_) : q1(q1_), q2(q2_) {}

tens_2q::tens_2q(std::complex<double> alpha00, std::complex<double> alpha01, 
                 std::complex<double> alpha10, std::complex<double> alpha11)
    : q1(alpha00, alpha10), q2(alpha01, alpha11) {}



// getters
qubit tens_2q::get_q1() const { return q1; }
qubit tens_2q::get_q2() const { return q2; }



// Affichage de l'état du système
void tens_2q::display() const {
    std::cout << "État du système à 2 qubits :" << std::endl;
    q1.display();
    q2.display();
}

// Produit scalaire
std::complex<double> tens_2q::operator|(const qubit_system& qs) const {
    const tens_2q* q = dynamic_cast<const tens_2q*>(&qs);
    if (q) {
        return (q1 | q->q1) * (q2 | q->q2);
    } else {
        throw std::runtime_error("Produit scalaire non défini pour ce type.");
    }
}

// Produit tensoriel
matrice tens_2q::operator&(const qubit_system& qs) const {
    const tens_2q* q = dynamic_cast<const tens_2q*>(&qs);
    if (q) {
        return (q1 & q->q1) * (q2 & q->q2);
    } else {
        throw std::runtime_error("Produit tensoriel non défini pour ce type.");
    }
}
