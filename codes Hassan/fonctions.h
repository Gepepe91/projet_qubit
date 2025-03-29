#ifndef FONCTIONS_H
#define FONCTIONS_H

#include <fstream>
#include <vector>
#include <complex>
#include <random>
#include "qubit.h"
#include "matrice.h"

using complexe = std::complex<double>;

qubit compute_derivative(qubit q, matrice H);
void rk4(qubit& q, matrice H, double dt);
void simulation_static(qubit q , matrice H , double dt , double T_init, double T_max, int n);
void simulation_dynamic(qubit q_init, double omega, double omega_0, double omega_1, double dt, double T_init, double T_max, int n);
qubit preparation_etat_plus(qubit q_init, double omega, double omega_0, double omega_1, double dt, double T_init, double T_max, int n);
void appliquer_bruit_aleatoire(qubit& q, double niveau_bruit, double* bruittheta, double* bruitphi);
void appliquer_bruit_aleatoire_ab(qubit& q, double niveau_bruit, double* bruitalpha, double* bruitbeta);
void simulation_avec_bruit(qubit q, matrice H, double dt, double T_init, double T_max, int n, double noise_level);
void simulation_oscillating_magnetic_field_with_noise_and_correction(qubit q_init, double omega, double omega_0, double omega_1, double dt, double T_init, double T_max, int n, double noise_level);
void simulation_multiple_noise_levels(qubit q_init, double omega, double omega_0, double omega_1,double dt, double T_init, double T_max, int n);
void simulation_multiple_noise_levels_repete(qubit q_init, double omega, double omega_0, double omega_1,
    double dt, double T_init, double T_max, int n, int repetitions);
complexe f(double t, complexe alpha, complexe beta, double omega, double omega_0, double omega_1);
complexe g(double t, complexe alpha, complexe beta, double omega, double omega_0, double omega_1);

#endif // FONCTIONS_H