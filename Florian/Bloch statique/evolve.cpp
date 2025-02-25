#include "evolve.h"
#include <cmath>
#include <iostream>

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


void Evolve::getCoordinates(const Qubit& qubit, double& x, double& y, double& z) {
    // Calcul des coordonnées à partir des amplitudes
    complex<double> psi0 = qubit.getPsi0();
    complex<double> psi1 = qubit.getPsi1();

    x = 2.0 * real(conj(psi0) * psi1);
    y = 2.0 * imag(conj(psi0) * psi1);
    z = qubit.prob_psi0() - qubit.prob_psi1();
}


void Evolve::runSimulation(const Qubit& qubit, double delta_t, double T_total, double gamma, double Bz, ofstream& outFile) {
    Qubit evolvingQubit = qubit ;  // Copie de l'état initial

    for (double t = 0; t < T_total; t += delta_t) {
        evolve_stat(evolvingQubit, gamma, Bz, delta_t);

        // Calcul des probabilités
        double prob0 = qubit.prob_psi0();
        double prob1 = qubit.prob_psi1();

        // Calcul des coordonnées de la sphère de Bloch
        double x, y, z;
        getCoordinates(evolvingQubit, x, y, z);

        // Écrire les résultats
    outFile << t << "," << evolvingQubit.prob_psi0() << "," << evolvingQubit.prob_psi1() << "," << x << "," << y << "," << z << ',' << acos(z/(sqrt(x*x+y*y+z*z))) <<"\n";
    }
}
