#ifndef CHAMP_H
#define CHAMP_H

#include <iostream>
#include <cmath>
#include <fstream>
#include "Matrice.h"

class H_Champ {
private:
    const double hbar = 1.05e-34;   // Constante de Planck réduite (J·s)
    const double q = -1.6e-19;      // Charge de l'électron (Coulombs)
    const double m = 9.11e-31;      // Masse de l'électron (kg)
    const double g = -2.0;          // Facteur gyromagnétique pour un électron
    double B;                       // Champ magnétique (Tesla)
    Matrice Sz;                     // Matrice de Pauli-Z
    Matrice Sx;                     // Matrice de Pauli-X
    Matrice Sy;                     // Matrice de Pauli-Y

    // Méthode générique pour calculer un Hamiltonien
    Matrice calculerHamiltonien(const Matrice& S, double coefficient) const;

public: 
    // Constructeur
    H_Champ(double B_) : B(B_), Sz(2), Sx(2), Sy(2) {
        Sz.PauliZ(); // Initialiser la matrice Pauli-Z
        Sx.PauliX(); // Initialiser la matrice Pauli-X
        Sy.PauliY(); // Initialiser la matrice Pauli-Y
    }

    double calculerGamma() const; // Méthode pour calculer γ, rapport gyromagnétique
    double calculerf() const;     // Méthode pour calculer fréquence

    Matrice calculerH() const;    // Méthode pour calculer le Hamiltonien H pour un champ statique
    Matrice calculerHt(double t,double omega) const;   // Méthode pour calculer le Hamiltonien H pour un champ variable

    // Méthode pour calculer Htot (la somme de l'Hamiltonien statique et dynamique) et enregistrer dans un fichier CSV
    void enregistrerHtot(double temps_max, double pas_temps,double omega) const;

};

#endif