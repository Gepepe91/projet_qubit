#include "evolve.h"
#include <cmath>
#include <iostream>
#include <complex>

using namespace std;

// Cette fonction effectue une évolution statique avec un champ Bz autour de l'axe Oz du qubit
void Evolve::evolve_stat(Qubit& qubit, double gamma, double Bz, double delta_t) {
    double omega_0 = - gamma * Bz;  // Fréquence de Larmor 

    // Matrice de rotation autour de l'axe Oz
    complex<double> U00 = exp(complex<double>(0, - omega_0 * delta_t / 2));
    complex<double> U01 = 0;
    complex<double> U10 = 0;
    complex<double> U11 = exp(complex<double>(0, omega_0 * delta_t / 2));

    // Application de la matrice de rotation
    complex<double> new_psi0 = U00 * qubit.getPsi0();
    complex<double> new_psi1 = U11 * qubit.getPsi1();

    qubit.setThetaPhi(2.0 * acos(abs(new_psi0)), arg(new_psi1 / new_psi0));
}

// Cette fonction calcule les coordonnées de la sphère de Bloch à partir des amplitudes du qubit
void Evolve::getCoordinates(const Qubit& qubit, double& x, double& y, double& z) {
    complex<double> psi0 = qubit.getPsi0();
    complex<double> psi1 = qubit.getPsi1();

    // Coordonnées sur la sphère de Bloch
    x = 2.0 * real(conj(psi0) * psi1);
    y = 2.0 * imag(conj(psi0) * psi1);
    z = qubit.prob_psi0() - qubit.prob_psi1();
}

// Méthode Runge-Kutta d'ordre 4 pour résoudre les équations différentielles couplées
void Evolve::runge_katta(Qubit& qubit, double gamma, double B1, double omega, double delta_t) {
    // Initialisation des amplitudes λ et μ
    complex<double> lambda = qubit.getPsi0();
    complex<double> mu = qubit.getPsi1();

    // Calcul de ω1
    double omega1 = - gamma * B1;

    // Calcul des dérivées à différents stades (k1, k2, k3, k4)
    complex<double> k1_lambda = -0.5 * omega1 * exp(complex<double>(0, (omega - gamma * B1) * 0)) * mu;
    complex<double> k1_mu = 0.5 * omega1 * exp(complex<double>(0, -(omega - gamma * B1) * 0)) * lambda;

    complex<double> k2_lambda = -0.5 * omega1 * exp(complex<double>(0, (omega - gamma * B1) * 0.5 * delta_t)) * (mu + 0.5 * delta_t * k1_mu);
    complex<double> k2_mu = 0.5 * omega1 * exp(complex<double>(0, -(omega - gamma * B1) * 0.5 * delta_t)) * (lambda + 0.5 * delta_t * k1_lambda);

    complex<double> k3_lambda = -0.5 * omega1 * exp(complex<double>(0, (omega - gamma * B1) * 0.5 * delta_t)) * (mu + 0.5 * delta_t * k2_mu);
    complex<double> k3_mu = 0.5 * omega1 * exp(complex<double>(0, -(omega - gamma * B1) * 0.5 * delta_t)) * (lambda + 0.5 * delta_t * k2_lambda);

    complex<double> k4_lambda = -0.5 * omega1 * exp(complex<double>(0, (omega - gamma * B1) * delta_t)) * (mu + delta_t * k3_mu);
    complex<double> k4_mu = 0.5 * omega1 * exp(complex<double>(0, -(omega - gamma * B1) * delta_t)) * (lambda + delta_t * k3_lambda);

    // Mise à jour des amplitudes λ et μ
    lambda += (delta_t / 6.0) * (k1_lambda + 2.0 * k2_lambda + 2.0 * k3_lambda + k4_lambda);
    mu += (delta_t / 6.0) * (k1_mu + 2.0 * k2_mu + 2.0 * k3_mu + k4_mu);

    // Mise à jour des états du qubit
    qubit.setThetaPhi(2.0 * acos(abs(lambda)), arg(mu / lambda));
}

// Méthode pour exécuter la simulation dynamique du qubit
void Evolve::runSimulation_dyn(Qubit& qubit, double gamma, double B1, double omega, double T_total, double delta_t, ofstream& outFile) {
    // Initialisation de l'état du qubit
    Qubit evolvingQubit = qubit;

    // Simulation temporelle
    for (double t = 0; t < T_total; t += delta_t) {
        // Applique Runge-Kutta pour évoluer l'état du qubit
        runge_katta(evolvingQubit, gamma, B1, omega, delta_t);

        // Calcul des probabilités et des coordonnées sur la sphère de Bloch
        double prob0 = evolvingQubit.prob_psi0();
        double prob1 = evolvingQubit.prob_psi1();
        double x, y, z;
        getCoordinates(evolvingQubit, x, y, z);

        // Écriture des résultats dans le fichier
        outFile << t << "," << prob0 << "," << prob1 << "," << x << "," << y << "," << z << "\n";
    }
}