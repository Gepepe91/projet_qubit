#ifndef EVOLVE_H
#define EVOLVE_H

#include "qubit.h"
#include <complex>
#include <fstream>

using namespace std;

struct Evolve {
    // Applique une évolution avec champ B statique // Oz au qubit (rotation autour de l'axe Oz)
    static void evolve_stat(Qubit& qubit, double gamma, double Bz, double delta_t);

    // Calcul des coordonnées sur la sphère de Bloch pour un qubit donné
    static void getCoordinates(const Qubit& qubit, double& x, double& y, double& z);

    // Exécute la simulation et écrit les résultats dans un fichier
    static void runSimulation(const Qubit& qubit, double delta_t, double T_total, double gamma, double Bz, ofstream& outFile);
    
    // Méthode Runge-Kutta d'ordre 4 pour résoudre les équations différentielles couplées
    static void runge_katta(Qubit& qubit, double gamma, double B1, double omega, double delta_t) ;

    // Méthode pour exécuter la simulation dynamique du qubit
    static void runSimulation_dyn(Qubit& qubit, double gamma, double B1, double omega, double T_total, double delta_t, ofstream& outFile) ;

};

#endif