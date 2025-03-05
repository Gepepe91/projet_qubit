#include <iostream>
#include <fstream>
#include <cmath>
#include "qubit.h"
#include "evolve.h"
#include "constant.h"

using namespace std ;

double degToRad(double degrees) {
    return degrees * M_PI / 180.0;
}

int main() {
    double Bz, theta_deg, phi_deg, delta_t, T_total;

    cout << "Entrez la valeur de Bz (en Tesla) : ";
    cin >> Bz;
    cout << "Entrez l'angle theta (en degres) : ";
    cin >> theta_deg;
    cout << "Entrez l'angle phi (en degres) : ";
    cin >> phi_deg;
    cout << "Entrez le temps total T_total (en secondes) : ";
    cin >> T_total;
    cout << "Entrez le pas de temps delta_t (en secondes) : ";
    cin >> delta_t;


    double theta = degToRad(theta_deg);
    double phi = degToRad(phi_deg);

    Qubit qubit(theta, phi);
    double gamma = Constants::g * Constants::q / (2 * Constants::m_e) ;
    double omega = - gamma * Bz ;

    ofstream outFile("results.csv");
    outFile << "t,prob_psi0,prob_psi1,x,y,z\n";

    Evolve::runSimulation(qubit, delta_t, T_total, gamma, Bz, outFile);

    outFile.close();
    cout << "Simulation terminee. Resultats enregistres dans results.csv\n";


    cout << "Periode de precession autour de OZ : " << 2*M_PI / omega << " s" << endl ;
    return 0;
}