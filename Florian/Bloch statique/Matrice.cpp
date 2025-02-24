#include "Matrice.h"
using namespace std; 

// Matrice associé à la porte Hadamard
void Matrice::gateHADAMARD(){
    matrice[0][0] = 1/sqrt(2) ;
    matrice[0][1] = 1/sqrt(2) ;
    matrice[1][0] = 1/sqrt(2) ;
    matrice[1][1] = -1/sqrt(2) ;
}


// Affichage de la matrice
void Matrice::affichermatrice() const {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                cout << matrice[i][j] << " ";
            }
            cout << endl;
        }
    }


// Matrice de Pauli suivant Oz
void Matrice::PauliZ() {
    matrice[0][0] = complex<double>(1,0) ;
    matrice[0][1] = complex<double>(0,0) ;
    matrice[1][0] = complex<double>(0,0) ;
    matrice[1][1] = complex<double>(-1,0) ;
}


// Matrice de Pauli suivant Ox
void Matrice::PauliX() {
    matrice[0][0] = complex<double>(0,0) ;
    matrice[0][1] = complex<double>(1,0) ;
    matrice[1][0] = complex<double>(1,0) ;
    matrice[1][1] = complex<double>(0,0) ;
}


// Matrice de Pauli suivant Ox
void Matrice::PauliY() {
    matrice[0][0] = complex<double>(0,0) ;
    matrice[0][1] = complex<double>(0,-1) ;
    matrice[1][0] = complex<double>(0,1) ;
    matrice[1][1] = complex<double>(0,0) ;
}

void Matrice::Identity() {
    matrice[0][0] = complex<double>(1,0) ;
    matrice[0][1] = complex<double>(0,0) ;
    matrice[1][0] = complex<double>(0,0) ;
    matrice[1][1] = complex<double>(1,0) ;
}



// Méthode pour modifier la valeur d'un élément de la matrice
void Matrice::setValue(int i, int j, complex<double> value) {
    matrice[i][j] = value;
}

// Méthode pour accéder à la valeur d'un élément de la matrice
complex<double> Matrice::getValue(int i, int j) const {
    return matrice[i][j];
}

// Méthode pour avoir la dimension de la matrice 
int Matrice::getDimension() const { 
    return dimension; 
    } 


// Surcharge pour l'addition 
Matrice Matrice::operator+(const Matrice& other) const {
    Matrice result(this->dimension); // Créer une nouvelle matrice pour stocker le résultat

    for (int i = 0; i < this->dimension; ++i) {
        for (int j = 0; j < this->dimension; ++j) {
            result.setValue(i, j, this->getValue(i, j) + other.getValue(i, j)); // Additionner les éléments correspondants
        }
    }

    return result;
}
