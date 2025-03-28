#include "qubit_system.h"
#include <iostream>

// Affichage par défaut (doit être redéfini par les classes filles)
void qubit_system::display() const {
    std::cout << "Affichage générique d'un système de qubits." << std::endl;
}

// Produit scalaire | (par défaut)
std::complex<double> qubit_system::operator|(const qubit_system& qs) const {
    throw std::runtime_error("Produit scalaire non défini pour ce type de système de qubits.");
}

// Produit tensoriel & (par défaut)
matrice qubit_system::operator&(const qubit_system& qs) const {
    throw std::runtime_error("Produit tensoriel non défini pour ce type de système de qubits.");
}
