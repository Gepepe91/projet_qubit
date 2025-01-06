#ifndef MATRICE_H
#define MATRICE_H

#include <vector>
#include <iostream>
#include <cmath>
#include <complex> 

using namespace std; 

class Matrice {

    private :
    vector<vector<complex<double>>> matrice ;
    int dimension ;


    public :

    const double hbar = 1.05e-34; // Constante de Planck réduite
    Matrice(int dimension_) : dimension(dimension_) {
        matrice.resize(dimension_, vector<complex<double>>(dimension_,complex<double>(0, 0))) ;
    } 

    void gateHADAMARD() ;
    void PauliZ() ;
    void PauliX() ;
    void PauliY() ;
    void Identity() ;
    void affichermatrice() const ;
    void setValue(int i, int j, complex<double> value);  
    complex<double> getValue(int i, int j) const;
    int getDimension() const ; 
    Matrice operator+(const Matrice& other) const;



};

#endif