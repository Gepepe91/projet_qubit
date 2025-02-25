#include<iostream>
#include"matrice.h"
#include<complex>
#include<cmath>
#include"qubit.h"

using complex = std::complex<double>; //ça améliore grandement la lisibilité

matrice::matrice() : mat(2 , std::vector<complex>(2 , complex(0,0))) {
    mat[0][0] = complex(1,0);
    mat[1][1] = complex(1,0);
    // 1 0
    // 0 1
}

matrice::matrice(complex a , complex b , complex c , complex d) : mat(2 , std::vector<complex>(2 , complex(0,0))) {
    mat[0][0] = a; mat[1][0] = c; mat[0][1] = b; mat[1][1] = d;
    // a b
    // c d
}

void matrice::set_as(complex a , complex b , complex c , complex d){
    mat[0][0] = a; mat[1][0] = c; mat[0][1] = b; mat[1][1] = d;
}

void matrice::display() {
    for ( std::vector<complex>& ligne : mat) {
        for ( complex& elem : ligne){
            std::cout << elem << "  ";
        }
        std::cout << std::endl;
    }
}

complex matrice::get_element(unsigned int ligne , unsigned int colonne) const {
    return mat[ligne][colonne];
}


void matrice::set_as_Pauli_x(){
    mat[0][0] = complex(0,0);
    mat[1][0] = complex(1,0);
    mat[0][1] = complex(1,0);
    mat[1][1] = complex(0,0);
}
void matrice::set_as_Pauli_y(){
    mat[0][0] = complex(0,0);
    mat[1][0] = complex(0,1);
    mat[0][1] = complex(0,-1);
    mat[1][1] = complex(0,0);
}
void matrice::set_as_Pauli_z(){
    mat[0][0] = complex(1,0);
    mat[1][0] = complex(0,0);
    mat[0][1] = complex(0,0);
    mat[1][1] = complex(-1,0);
}

matrice matrice::operator*(const complex z){
    return matrice(
        mat[0][0] * z, mat[0][1] * z,
        mat[1][0] * z, mat[1][1] * z
    );
}