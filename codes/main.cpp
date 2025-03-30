#include<iostream>
#include<cmath>
#include<complex>
#include<random>

#include"matrice.h"
#include"qubit_system.h"
#include"qubit.h"
#include"tens_2q.h"
#include"tens_Nq.h"

#include"constants.h"
#include"include/fonctions.h"

using complexe = std::complex<double>;
using namespace std;


int main(){


// définition des grandeurs physiques :

double B_z = 0.01 ; //champ B0 de 10 mT
double B_xy = 0.001;

double omega_0 = - gamma_e * B_z ; //fréquence de Larmor
double omega_1 = - gamma_e*B_xy;
omega_1 = std::abs(omega_1); // Ensure omega_1 is positive
//double omega = 2.288e9; //tel que w - w0 = 3*w1 et amplitude varie entre 0.9-1 et 0-0.1
double omega = omega_0; //valeur de w à résonance

// grandeurs de la simulation :

int n = 1000;
double T_1 = (2*M_PI)/omega_1;
double T_init = 0. ;
double T_max = T_init + 5*T_1; // arbitraire aussi, temps maximale de la simulation
double dt = (T_max-T_init) / n; //incrément de temps arbitraire comme étant une fraction de la période
double noise_level = 0.05; // Niveau de bruit

// Hamiltonien de précession de phase

matrice H_null(0,0,0,0);
matrice H_stat_Oz(hbar*omega_0/2. , 0 , 0. , -hbar*omega_0/2.); //Champ selon Oz uniquement
matrice H_stat_Ox_Oz(hbar*omega_0/2. , hbar*omega_1/2 , hbar*omega_1/2 ,-hbar*omega_0/2. ); //Champ selon 0x et Oz
//matrice H_stat_Ox_Oz_coupled(omega_0. , omega_1 , omega_1 ,omega_0 );

////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                        /*Utilisation de la classe qubit*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////

qubit q;                                                        //Constructeur par défaut
qubit q1(0,0);                                                  //Initialisation avec theta et phi
qubit q2(complexe (1,0) , complexe (1./sqrt(2) , 1./sqrt(2)));  //Initialisation avec alpha et beta

////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                        /*Oscillations de Rabi non bruitées*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////

simulation_static(q1 , H_stat_Oz , dt , T_max); //simulation de l'évolution d'un qubit dans un champs constant
simulation_dynamic(q1 , omega , omega_0 , omega_1 , dt , T_init , T_max , n); // simulation de l'évolution d'un qubit dans un champ dynamique

////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                        /*Préparation de qubit simples dans l'état |+>*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////

preparation_etat_plus(q1 , omega , omega_0 , omega_1 , dt , T_init , T_max, n);
preparation_bruit_correction(q1 , omega , omega_0 , omega_1 , dt , T_init , T_max, n , noise_level);

////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                        /*Calculs de l'information de Fisher>*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////

//fisher_1qubit();
//fisher_2qubit_intrique();


    return 0;
}