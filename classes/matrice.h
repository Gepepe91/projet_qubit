#ifndef MATRICE_H
#define MATRICE_H

#include <vector>
#include <iostream>
#include <stdexcept>
#include <complex>
#include "qubit.h"

using complexe = std::complex<double>;

class matrice {
private:
    std::vector<std::vector<complexe>> mat; // Représentation de la matrice sous forme de vecteur de vecteurs
    size_t lignes; // Nombre de lignes
    size_t colonnes; // Nombre de colonnes

public:
    // Constructeurs
    matrice(); // Constructeur par défaut, matrice identité 2x2
    matrice(size_t lignes, size_t colonnes);
    matrice(complexe a, complexe b, complexe c, complexe d); // Constructeur pour matrice 2x2 avec 4 éléments
    
    // Modificateurs
    void set_as(complexe a, complexe b, complexe c, complexe d); // Modifie la matrice avec 4 valeurs
    void set_as_Pauli_x(); // Matrice Pauli X
    void set_as_Pauli_y(); // Matrice Pauli Y
    void set_as_Pauli_z(); // Matrice Pauli Z
    void setElement(size_t i, size_t j, double valeur); // Modifie un élément de la matrice
    
    // Getters
    complexe get_element(unsigned int ligne, unsigned int colonne) const; // Accède à un élément de la matrice
    
    // Affichage
    void display() const; // Affiche la matrice

    // Opérations sur matrices
    matrice operator*(const complexe z) const; // Produit de la matrice avec un scalaire
    matrice operator*(const matrice& m) const; // Produit de matrices
    qubit operator*(const qubit& q) const; // Produit de la matrice avec un qubit

    matrice operator+(const complexe z) const; // Addition de la matrice avec un scalaire
    matrice operator+(const matrice& n) const; // Addition de matrices

    matrice dephasage(double xi); // Matrice de déphasage
};

#endif // MATRICE_H
