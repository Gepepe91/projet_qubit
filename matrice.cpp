#include<iostream>
#include"matrice.h"
#include"complexe.h"
#include<cmath>

matrice::matrice() : mat(2 , std::vector<complexe>(2 , complexe(0,0))) {
    mat[0][0] = complexe(1,0);
    mat[1][1] = complexe(1,0);
    // 1 0
    // 0 1
}

matrice::matrice(complexe a , complexe b , complexe c , complexe d) : mat(2 , std::vector<complexe>(2 , complexe(0,0))) {
    mat[0][0] = a; mat[1][0] = b; mat[0][1] = c; mat[1][1] = d;
    // a c
    // b d
}

void matrice::display() {
    for ( std::vector<complexe>& ligne : mat) {
        for ( complexe& elem : ligne){
            elem.afficher(); std::cout << "  ";
        }
        std::cout << std::endl;
    }
}

complexe matrice::get_element(unsigned int ligne , unsigned int colonne) const {
    return mat[ligne][colonne];
}


void matrice::set_as_Hadamard(){
    mat[0][0] = complexe(1/sqrt(2),0);
    mat[1][0] = complexe(1/sqrt(2),0);
    mat[0][1] = complexe(1/sqrt(2),0);
    mat[1][1] = complexe(-1/sqrt(2),0);
}
void matrice::set_as_Pauli_x(){
    mat[0][0] = complexe(0,0);
    mat[1][0] = complexe(1,0);
    mat[0][1] = complexe(1,0);
    mat[1][1] = complexe(0,0);
}
void matrice::set_as_Pauli_y(){
    mat[0][0] = complexe(0,0);
    mat[1][0] = complexe(0,1);
    mat[0][1] = complexe(0,-1);
    mat[1][1] = complexe(0,0);
}
void matrice::set_as_Pauli_z(){
    mat[0][0] = complexe(1,0);
    mat[1][0] = complexe(0,0);
    mat[0][1] = complexe(0,0);
    mat[1][1] = complexe(-1,0);
}