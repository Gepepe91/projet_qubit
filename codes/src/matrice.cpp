#include<iostream>
#include "matrice.h"
#include "qubit.h"

matrice::matrice() : mat(2, std::vector<complexe>(2, complexe(0,0))) {
    mat[0][0] = complexe(1,0);
    mat[1][1] = complexe(1,0);
}

matrice::matrice(complexe a, complexe b, complexe c, complexe d) : mat(2, std::vector<complexe>(2, complexe(0,0))) {
    mat[0][0] = a; mat[0][1] = b;
    mat[1][0] = c; mat[1][1] = d;
}

matrice::matrice(size_t lignes, size_t colonnes) : lignes(lignes), colonnes(colonnes) {
    mat.resize(lignes, std::vector<complexe>(colonnes, complexe(0, 0)));
}

void matrice::set_as(complexe a, complexe b, complexe c, complexe d) {
    mat[0][0] = a; mat[0][1] = b;
    mat[1][0] = c; mat[1][1] = d;
}

void matrice::setElement(size_t i, size_t j, double valeur) {
    if (i >= lignes || j >= colonnes) {
        throw std::out_of_range("Indices out of bounds.");
    }
    mat[i][j] = valeur;
}

void matrice::display() const {
    for (const auto& ligne : mat) {
        for (const auto& elem : ligne) {
            std::cout << elem << "  ";
        }
        std::cout << std::endl;
}}

complexe matrice::get_element(unsigned int ligne, unsigned int colonne) const {
    return mat[ligne][colonne];
}

void matrice::set_as_Pauli_x() {
    set_as(0, 1, 1, 0);
}

void matrice::set_as_Pauli_y() {
    set_as(0, complexe(0,-1), complexe(0,1), 0);
}

void matrice::set_as_Pauli_z() {
    set_as(1, 0, 0, -1);
}

matrice matrice::operator*(const complexe z) const {
    return matrice(mat[0][0] * z, mat[0][1] * z,
                   mat[1][0] * z, mat[1][1] * z);
}

matrice matrice::operator*(const matrice& m) const {
    return matrice(
        mat[0][0] * m.mat[0][0] + mat[0][1] * m.mat[1][0],
        mat[0][0] * m.mat[0][1] + mat[0][1] * m.mat[1][1],
        mat[1][0] * m.mat[0][0] + mat[1][1] * m.mat[1][0],
        mat[1][0] * m.mat[0][1] + mat[1][1] * m.mat[1][1]
    );
}

qubit matrice::operator*(const qubit& q) const {
    return qubit(
        mat[0][0] * q.get_alpha() + mat[0][1] * q.get_beta(),
        mat[1][0] * q.get_alpha() + mat[1][1] * q.get_beta()
    );
}

matrice matrice::operator+(const complexe z) const {
    return matrice(mat[0][0] + z, mat[0][1] + z,
                   mat[1][0] + z, mat[1][1] + z);
}

matrice matrice::operator+(const matrice& n) const {
    return matrice(mat[0][0] + n.mat[0][0], mat[0][1] + n.mat[0][1],
                   mat[1][0] + n.mat[1][0], mat[1][1] + n.mat[1][1]);
}

matrice matrice::dephasage(double xi) {
    return matrice(
        std::cos(xi), -std::sin(xi),
        std::sin(xi), std::cos(xi)
    );
}
