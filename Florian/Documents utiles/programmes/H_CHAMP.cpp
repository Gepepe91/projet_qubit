#include "H_CHAMP.h"
#include <fstream>
#include <complex>

using namespace std ;

// Méthode pour calculer γ
double H_Champ::calculerGamma() const {
    return g * q / (2 * m);
}

// Méthode pour calculer la fréquence
double H_Champ::calculerf() const {
    return calculerGamma() * B;
}

// Méthode générique pour calculer un Hamiltonien
Matrice H_Champ::calculerHamiltonien(const Matrice& S, double coefficient) const {
    Matrice H(S.getDimension()); // Taille identique à Sz
    for (int i = 0; i < H.getDimension(); ++i) {
        for (int j = 0; j < H.getDimension(); ++j) {
            H.setValue(i, j, coefficient * S.getValue(i, j));
        }
    }
    return H;
}

// Méthode pour calculer le Hamiltonien H statique
Matrice H_Champ::calculerH() const {
    double w0 = calculerf(); // Pulsation cyclotron
    return calculerHamiltonien(Sz, 0.5 * hbar * w0);
}

// Méthode pour calculer le Hamiltonien H variable à un instant t
Matrice H_Champ::calculerHt(double t,double omega) const {
    double w1 = calculerf(); // Fréquence de Rabi
    return calculerHamiltonien(Sx, -0.5 * hbar * w1 * cos(omega * t)) +  
           calculerHamiltonien(Sy, 0.5 * hbar * w1 * sin(omega * t));   
}

// Méthode pour calculer Htot (la somme de l'Hamiltonien statique et dynamique) et enregistrer dans un fichier CSV
void H_Champ::enregistrerHtot(double temps_max, double pas_temps, double omega) const {
    // Ouvre un fichier pour écrire les résultats
    ofstream fichier("Htot.csv");
    
    // En-tête du fichier CSV
    fichier << "Temps (s), H_11_real, H_11_imag, H_12_real, H_12_imag, H_21_real, H_21_imag, H_22_real, H_22_imag\n\n";  
    
    // Calcul de Htot pour chaque valeur de t de 0 à temps_max avec un pas de pas_temps
    for (double t = 0; t <= temps_max; t += pas_temps) {
        Matrice Htot = calculerH() + calculerHt(t,omega);  // Somme des Hamiltoniens statique et dynamique à l'instant t
        
        // Récupérer les valeurs de la matrice Htot
        complex<double> H_11 = Htot.getValue(0, 0);
        complex<double> H_12 = Htot.getValue(0, 1);
        complex<double> H_21 = Htot.getValue(1, 0);
        complex<double> H_22 = Htot.getValue(1, 1);
        
        // Enregistrer les parties réelle et imaginaire de chaque élément de la matrice
        fichier << t << ", "
                << H_11.real() << ", " << H_11.imag() << ", "  // Partie réelle et imaginaire
                << H_12.real() << ", " << H_12.imag() << ", "
                << H_21.real() << ", " << H_21.imag() << ", "
                << H_22.real() << ", " << H_22.imag() << "\n";
    }
    
    fichier.close();  // Fermer le fichier après l'écriture
}