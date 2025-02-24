#include <iostream>
#include "H_CHAMP.h"

using namespace std;

int main() {

    // Initialisation des champs magnétiques statique et variable
    double B_static = 5.0; // Champ magnétique statique en Tesla
    double B_variable = 2.0; // Champ magnétique variable en Tesla


    // Création d'objets H_Champ pour les deux configurations
    H_Champ h_static(B_static);
    H_Champ h_variable(B_variable);
    double omega = 1.0; // Fréquence du champ B périodique en rad/s (exemple)


    // Calcul des valeurs physiques et affichage
    cout << "Gamma (rapport gyromagnetique) : " << h_static.calculerGamma() << " 1/s.T" << endl;
    cout << "Frequence de Larmor (pour B = " << B_static << " T) : " << h_static.calculerf() << " Hz" << endl;
    cout << "Frequence de Rabi (pour B = " << B_variable << " T) : " << h_variable.calculerf() << " Hz" << endl;


    // Enregistrer la somme des Hamiltoniens (Htot) dans un fichier CSV
    double temps_max = 10.0; // Temps maximum (en secondes)
    double pas_temps = 0.1;  // Pas de temps (en secondes)
    h_variable.enregistrerHtot(temps_max, pas_temps,omega); // Enregistrement des données dans le fichier "Htot.csv"

    cout << "\nLes donnees de Htot ont ete enregistrees dans 'Htot.csv'" << endl;

    return 0;
}