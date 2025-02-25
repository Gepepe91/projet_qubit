#ifndef MATRICE_H
#define MATRICE_H

/*
Je m'inspire de la classe Matrice de Florian, mais en essayant de généraliser un peu

Intêret de créer un objet matrice :
- permet de définir toutes les matrices qui nous intéresse : Pauli, Hadamard
Une façon originale : méthodes .init_as_Hadamard() , .init_as_Pauli_x, .init_as_Pauli_y
- méthode privée qui normalise
- dans la classe qubit, il faut ajouter une méthode qubit.transformation(matrice m) qui transforme le qubit selon la matrice
--> en somme, la classe matrice est assez simple, elle n'est pas directement responsable des transformations
Il est donc possible de s'inspirer très fortement du cours

D'ailleurs, c'est plus une classe d'isomorphisme qu'une classe de matrice : que des matrices 2x2 !
*/

#include<complex>
#include<vector>
#include<iostream>
#include<cmath>
#include"qubit.h"

class qubit;

using complex = std::complex<double>; //ça améliore grandement la lisibilité

class matrice{

    private:

    std::vector<std::vector<complex>> mat; // la matrice de Pauli y a des coef complexes

    public:

    matrice(); // le constructeur par défaut doit donner l'identité 2x2;
    matrice(complex a , complex c , complex b , complex d);
    // a b
    // c d

    void display();
    complex get_element(unsigned int ligne , unsigned int colonne) const;

    void set_as(complex a , complex b , complex c , complex d);

    void set_as_Pauli_x();
    void set_as_Pauli_y();
    void set_as_Pauli_z();

    //Surcharge opérateurs

    matrice operator*(const complex z); // matrice * complexe
    qubit operator*(const qubit q);     // matrice * qubit


};

#endif