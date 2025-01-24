#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include "qubit.h"
#include "matrice.h"

double DT=0.01;   // Pas de temps
double T_MAX=30;  // Temps maximal de la simulation

using complexe = std::complex<double>;
using namespace std;

// Fonction pour calculer f(t, psi) = -i * H * psi , ici correspond à la dérivée (éq de schrödinger)
// cette fonction revoie un nouveau qubit !
qubit compute_derivative(qubit q, matrice H) {
    // Appliquer la matrice H au qubit
    qubit new_q = q;
    new_q.transform(H);
    // Calculer la dérivée de l'état en utilisant l'hamiltonien (multiplication par -i/hbar)
    complexe temp_alpha(-new_q.get_alpha().imag() , new_q.get_alpha().real()) ; // -i * alpha
    complexe temp_beta(-new_q.get_beta().imag() , new_q.get_beta().real()) ; // -i *beta
    new_q.set_alpha(temp_alpha);
    new_q.set_beta(temp_beta);
    return new_q;
}

// Implémentation de la méthode Runge-Kutta de 4e ordre
void rk4(qubit& q, matrice H, double dt) {
    
    complexe a , b ; //le constructeur des complexes par défaut donne 0 + Oi

    //initialisation à alpha = 0 , beta = 0 des ki de RK4
    qubit k1(a , b);
    qubit k2(a , b);
    qubit k3(a , b);
    qubit k4(a , b);

    // k1
    k1 = compute_derivative(q, H);
    // k2
    k2 = compute_derivative(q + k1*(0.5*dt), H);
    // k3
    k3 = compute_derivative(q + k2*(0.5*dt) , H);
    // k4
    k4 = compute_derivative(q + k3*dt , H);

    // Mise à jour de l'état du qubit
    q = q + (k1 + k2*2 + k3*2 + k4) * (dt/6); // l'ordre est un peu chiant, difficulté de la surcharge en C++

}

int main() {
    // Hamiltonien de précession de phase (sans perturbation)
    double omega_0 = 1.0;
    matrice H(omega_0/2. , 0 , 0. , -omega_0/2.);

    H.display();

    // Initialisation de l'état |psi(0)> = |0>
    complexe alpha0(1,0) , beta0(1,0);
    qubit q(alpha0 , beta0);

    //Ouverture du fichier
    ofstream fichier;
    fichier.open("data.csv");

    // Simulation
    fichier << "Time\talpha\tbeta" << endl;//Re(alpha)\tIm(alpha)\tRe(beta)\tIm(beta)" << endl;
    for (double t = 0; t <= T_MAX; t += DT) {
        fichier << t << "\t" << q.get_alpha() << "\t" << q.get_beta() << endl;// << "\t" << q.get_alpha().real() << "\t" << q.get_alpha().imag() << "\t"
                  //<< q.get_beta().real() << "\t" << q.get_beta().imag() << std::endl;
        rk4(q, H, DT); // Mise à jour de |psi> = notre qubit avec RK4
    }

    //Fermeture du fichier
    fichier.close();

    return 0;
}
