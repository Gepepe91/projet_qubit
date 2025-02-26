#include<iostream>
#include<cmath>
#include"matrice.h"
#include<complex>
#include"qubit.h"
#include"constants.h"
#include"fonctions.h"

using complexe = std::complex<double>;
using namespace std;



int main(){

    // définition des grandeurs physiques :

    double B_z = 0.01 ; //champ B0 de 10 mT
    double B_xy = 0.01;

    double omega_0 = - gamma_e * B_z ; //fréquence de Larmor
    double omega_1 = - gamma_e*B_xy;

    double omega = omega_0; //arbitraire, controlé par l'expérimentateur, résonnance etc

    // grandeurs de la simulation :

    double T_0 = (2*M_PI)/omega_0;
    double dt = T_0 / 100000; //incrément de temps arbitraire comme étant une fraction de la période
    double T_max = 10*T_0; // arbitraire aussi, temps maximale de la simulation

    // Hamiltonien de précession de phase

    matrice H_null(0,0,0,0);
    matrice H_stat_Oz(hbar*omega_0/2. , 0 , 0. , -hbar*omega_0/2.); //Champ selon Oz uniquement
    matrice H_stat_Ox_Oz(hbar*omega_0/2. , hbar*omega_1/2 , hbar*omega_1/2 ,-hbar*omega_0/2. ); //Champ selon 0x et Oz

    // Initialisation de notre qubit , état initial

    qubit q(-2*M_PI,0); // theta et phi

    //simulation

    simulation(q , H_stat_Ox_Oz , dt , T_max);

    

    return 0;
}