#include<iostream>
#include<cmath>
#include"matrice.h"
#include<complex>
#include <random>
#include"qubit.h"
#include"constants.h"
#include"fonctions.h"

using complexe = std::complex<double>;
using namespace std;



int main(){

/*    // définition des grandeurs physiques :

    double B_z = 0.01 ; //champ B0 de 10 mT
    double B_xy = 0.001;

    double omega_0 = - gamma_e * B_z ; //fréquence de Larmor
    double omega_1 = - gamma_e*B_xy;

    //double omega = 2.288e9; //tel que w - w0 = 3*w1 et amplitude varie entre 0.9-1 et 0-0.1
    double omega = omega_0; //valeur de w à résonance

    // grandeurs de la simulation :

    int n = 1000000;
    double T_1 = (2*M_PI)/omega_1;
    double T_init = 0. ;
    double T_max = T_init + 3*T_1; // arbitraire aussi, temps maximale de la simulation
    double dt = (T_max-T_init) / n; //incrément de temps arbitraire comme étant une fraction de la période

    // Hamiltonien de précession de phase

    matrice H_null(0,0,0,0);
    matrice H_stat_Oz(hbar*omega_0/2. , 0 , 0. , -hbar*omega_0/2.); //Champ selon Oz uniquement
    matrice H_stat_Ox_Oz(hbar*omega_0/2. , hbar*omega_1/2 , hbar*omega_1/2 ,-hbar*omega_0/2. ); //Champ selon 0x et Oz
    //matrice H_stat_Ox_Oz_coupled(omega_0. , omega_1 , omega_1 ,omega_0 );

    // Initialisation de notre qubit , état initial

    qubit q(0,0); // theta et phi

    //simulation

    simulation_static(q , H_stat_Oz , dt , T_max);

    simulation_dynamic(q  , omega , omega_0 , omega_1 ,  dt ,T_init, T_max , n);

    //Maintenant, on veut faire une fonction qui permet de préparer notre qubit dans l'état |+> = (|0> + |1>)/sqrt(2).
    //En pratique, cette fonction ressemble beaucoup à simulation_dynamic(), sauf que on stoppe l'itération lorsque on
    //est très proche de norm(alpha) = norm(beta) = 0.5.
    // !!! preparation_plus() renvoie un qubit, ce n'est pas une fonction void.

    qubit q_plus;

    q_plus = preparation_etat_plus(q  , omega , omega_0 , omega_1 ,  dt ,T_init, T_max , n);

    cout << q_plus.get_abs_alpha2() << endl;
    cout << q_plus.get_abs_beta2() << endl;

    */

    one_qubit(200 , M_PI/4.);

    two_qubits(100 ,  M_PI/4.);
    

    return 0;
}