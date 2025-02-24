#include <iostream>
#include <fstream>
#include <cmath>
#include "qubit.h"
#include "evolve.h"
#include "constant.h"
#include <string>


using namespace std;

double degToRad(double degrees) {
    return degrees * M_PI / 180.0;
}

int main() {
    double Bz, B1, omega, theta_deg, phi_deg, delta_t, T_total;
    string nom ;

    // Demander les entrées à l'utilisateur
    cout << "Entrez la valeur de Bz (en Tesla) : ";
    cin >> Bz;
    cout << "Entrez la valeur de B1 (en Tesla) : ";
    cin >> B1;
    cout << "Entrez la fréquence d'oscillation w (en Hz) : ";
    cin >> omega;
    cout << "Entrez l'angle theta (en degres) : ";
    cin >> theta_deg;
    cout << "Entrez l'angle phi (en degres) : ";
    cin >> phi_deg;
    cout << "Entrez le temps total T_total (en secondes) : ";
    cin >> T_total;
    cout << "Entrez le pas de temps delta_t (en secondes) : ";
    cin >> delta_t;

    // Conversion des angles en radians
    double theta = degToRad(theta_deg);
    double phi = degToRad(phi_deg);

    // Initialisation du qubit avec les paramètres theta et phi
    Qubit qubit(theta, phi);

    // Calcul du gamma
    double gamma = - Constants::g * Constants::q / (2 * Constants::m_e); 

    // Ouvrir le fichier pour écrire les résultats
    cout << "Nom du fichier .csv : " ;
    cin >> nom ;
    ofstream outFile(nom);
    outFile << "t,prob_psi0,prob_psi1,x,y,z\n";

    // Appel de la méthode d'évolution dynamique en utilisant Runge-Kutta
    Evolve::runSimulation_dyn(qubit, gamma, B1, omega, T_total, delta_t, outFile);

    // Fermer le fichier de résultats
    outFile.close();
    cout << "Simulation terminee. Resultats enregistres\n";

    return 0;
}