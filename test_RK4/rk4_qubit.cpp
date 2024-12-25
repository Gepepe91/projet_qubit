#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include "qubit.h"
#include "matrice.h"

double DT=0.01;   // Pas de temps
double T_MAX=10;  // Temps maximal de la simulation

using complexe = std::complex<double>;
using namespace std;

// Fonction pour calculer f(t, psi) = -i * H * psi
void compute_derivative(qubit q, matrice H, qubit dpsi) {
    // Appliquer la matrice H au qubit
    q.transform(H);

    // Calculer la dérivée de l'état en utilisant l'hamiltonien (multiplication par -i)
    complexe dpsi_alpha(-q.get_alpha().imag() , q.get_alpha().real()) ; // -i * alpha
    complexe dpsi_beta(-q.get_beta().imag() , q.get_beta().real()) ; // -i *beta
    dpsi.set_alpha(dpsi_alpha) ; dpsi.set_beta(dpsi_beta);
    //qubit dpsi(dpsi_alpha , dpsi_beta); //construction de dpsi
}

// Implémentation de la méthode Runge-Kutta de 4e ordre
void runge_kutta4(qubit &q, matrice H, double dt) {
    
    complexe a , b ; //le constructeur des complexes par défaut donne 0 + Oi

    qubit k1(a , b);
    qubit k2(a , b);
    qubit k3(a , b); //initialisation à alpha = 0 , beta = 0
    qubit k4(a , b);

    // k1
    compute_derivative(q, H, k1);

    // k2
    complexe k2_q_temp_alpha(q.get_alpha().real() + 0.5 * dt * k1.get_alpha().real() , q.get_alpha().imag() + 0.5 * dt * k1.get_alpha().imag());
    complexe k2_q_temp_beta(q.get_beta().real() + 0.5 * dt * k1.get_beta().real() , q.get_beta().imag() + 0.5 * dt * k1.get_beta().imag());
    
    qubit k2_q_temp(k2_q_temp_alpha , k2_q_temp_beta);

    compute_derivative(k2_q_temp, H, k2);

    // k3
    complexe k3_q_temp_alpha(q.get_alpha().real() + 0.5 * dt * k1.get_alpha().real() , q.get_alpha().imag() + 0.5 * dt * k1.get_alpha().imag());
    complexe k3_q_temp_beta(q.get_beta().real() + 0.5 * dt * k1.get_beta().real() , q.get_beta().imag() + 0.5 * dt * k1.get_beta().imag());
    
    qubit k3_q_temp(k3_q_temp_alpha , k3_q_temp_beta);

    compute_derivative(k3_q_temp , H , k3);

    // k4
    complexe k4_q_temp_alpha(q.get_alpha().real() + 0.5 * dt * k1.get_alpha().real() , q.get_alpha().imag() + 0.5 * dt * k1.get_alpha().imag());
    complexe k4_q_temp_beta(q.get_beta().real() + 0.5 * dt * k1.get_beta().real() , q.get_beta().imag() + 0.5 * dt * k1.get_beta().imag());
    
    qubit k4_q_temp(k4_q_temp_alpha , k4_q_temp_beta);

    compute_derivative(k4_q_temp , H , k4);

    // Mise à jour de l'état du qubit q
    complexe q_alpha(q.get_alpha().real() + (dt / 6.0) * (k1.get_alpha().real() + 2.0 * k2.get_alpha().real() + 2.0 * k3.get_alpha().real() + k4.get_alpha().real()),
                     q.get_alpha().imag() + (dt / 6.0) * (k1.get_alpha().imag() + 2.0 * k2.get_alpha().imag() + 2.0 * k3.get_alpha().imag() + k4.get_alpha().imag()));
    complexe q_beta(q.get_beta().real() + (dt / 6.0) * (k1.get_beta().real() + 2.0 * k2.get_beta().real() + 2.0 * k3.get_beta().real() + k4.get_beta().real()),
                     q.get_beta().imag() + (dt / 6.0) * (k1.get_beta().imag() + 2.0 * k2.get_beta().imag() + 2.0 * k3.get_beta().imag() + k4.get_beta().imag()));
    q = qubit(q_alpha,q_beta);

}

int main() {
    // Hamiltonien de précession de phase (sans perturbation)
    double omega_0 = 1.0;
    matrice H(omega_0/2. , 0. , 0. , -omega_0 /2.);

    // Initialisation de l'état |psi(0)> = |0>
    complexe alpha0(1,0) , beta0(0,0);
    qubit q(alpha0 , beta0);

    //Ouverture du fichier
    ofstream fichier;
    fichier.open("data.csv");

    // Simulation
    fichier << "Time\tRe(psi0)\tIm(psi0)\tRe(psi1)\tIm(psi1)" << endl;
    for (double t = 0; t <= T_MAX; t += DT) {
        fichier << t << "\t" << q.get_alpha().real() << "\t" << q.get_alpha().imag() << "\t"
                  << q.get_beta().real() << "\t" << q.get_beta().imag() << std::endl;
        runge_kutta4(q, H, DT); // Mise à jour de |psi> avec RK4
    }

    //Fermeture du fichier
    fichier.close();

    return 0;
}
