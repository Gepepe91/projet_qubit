#include<iostream>
#include<cmath>
#include<fstream>
#include<complex>
#include"qubit.h"
#include"matrice.h"

using complexe = std::complex<double>;
using namespace std;


double DT=0.05;   // Pas de temps
double T_MAX=10;  // Temps maximal de la simulation

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



int main(){

    // Hamiltonien de précession de phase (sans perturbation)
    double omega_0 = 1.0;
    matrice H(omega_0/2. , 0 , 0. , -omega_0/2.);

    // Initialisation de notre qubit , état initial
    complexe alpha0(sqrt(2),sqrt(2)) , beta0(sqrt(2),sqrt(2));
    qubit q(alpha0 , beta0);

    //Ouverture du fichier
    ofstream fichier;
    fichier.open("data.csv");

    // Simulation
    fichier << "Time,theta,phi" << endl;
    for (double t = 0; t <= T_MAX; t += DT) {
        fichier << t  << "," << q.get_theta() << "," << q.get_phi() << endl;
        rk4(q, H, DT); // Mise à jour de |psi> = notre qubit avec RK4
        q.normalize(); // renormalisation, sinon le qubit se dénormalise au cours du temps (erreur numérique)
    }

    //Fermeture du fichier
    fichier.close();

    return 0;
}
