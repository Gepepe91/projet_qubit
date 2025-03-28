#ifndef MATRICE_H
#define MATRICE_H

#include <vector>
#include <iostream>
#include <stdexcept>
#include <complex>
//#include "qubit.h" Ca foutait la merdre
class qubit; //déclaration anticipée

using complexe = std::complex<double>;

class matrice {
private:
    std::vector<std::vector<complexe>> mat; // Matrix represented as a vector of vectors
    size_t lignes; // Number of rows
    size_t colonnes; // Number of columns

public:
    // Constructors
    matrice(); // Default constructor, identity matrix 2x2
    matrice(size_t lignes, size_t colonnes);
    matrice(complexe a, complexe b, complexe c, complexe d); // 2x2 matrix with 4 elements
    
    // Modifiers
    void set_as(complexe a, complexe b, complexe c, complexe d); // Modifies matrix with 4 values
    void set_as_Pauli_x(); // Pauli X matrix
    void set_as_Pauli_y(); // Pauli Y matrix
    void set_as_Pauli_z(); // Pauli Z matrix
    void setElement(size_t i, size_t j, double valeur); // Modify matrix element
    
    // Getters
    complexe get_element(unsigned int ligne, unsigned int colonne) const; // Access matrix element
    
    // Display
    void display() const; // Display matrix

    // Matrix operations
    matrice operator*(const complexe z) const; // Matrix multiplication with a scalar
    matrice operator*(const matrice& m) const; // Matrix multiplication
    qubit operator*(const qubit& q) const; // Matrix multiplication with a qubit

    matrice operator+(const complexe z) const; // Matrix addition with a scalar
    matrice operator+(const matrice& n) const; // Matrix addition with another matrix

    matrice dephasage(double xi); // Dephasing matrix
};

#endif
