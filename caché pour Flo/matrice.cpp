#include<iostream>
#include"matrice.h"
#include<complex>
#include<cmath>

using complex = std::complex<double>; //ça améliore grandement la lisibilité

matrice::matrice() : mat(2 , std::vector<complex>(2 , complex(0,0))) {
    mat[0][0] = complex(1,0);
    mat[1][1] = complex(1,0);
    // 1 0
    // 0 1
}

matrice::matrice(complex a , complex b , complex c , complex d) : mat(2 , std::vector<complex>(2 , complex(0,0))) {
    mat[0][0] = a; mat[1][0] = b; mat[0][1] = c; mat[1][1] = d;
    // a c
    // b d
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


void matrice::set_as_Hadamard(){
    mat[0][0] = complex(1/sqrt(2),0);
    mat[1][0] = complex(1/sqrt(2),0);
    mat[0][1] = complex(1/sqrt(2),0);
    mat[1][1] = complex(-1/sqrt(2),0);
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