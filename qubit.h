#ifndef QUBIT_H
#define QUBIT_H

#include"complexe.h"

class qubit

/*
un qubit peut se modéliser de cette façon :
|psy(t)> = alpha|0> + beta|1> , alpha et beta complexes --> utilisation de la classe complexes créé en TD

avec |alpha|**2 + |beta|**2 = 1

un qubit est un objet relativement abstrait

1ere classe très générale, attributs non physiques
--> classes filles étant des qubits "réels"

les deux attributs d'un qubit sont alpha et beta, respectivement associés à |0> et |1>

il y a aussi deux autres attributs utiles pour la représentation selon la sphère de Bloch : theta et phi
--> méthodes privées de synchronisation entre ces deux choses représentations

*/

 {

private:

complexe alpha , beta;
double theta , phi;

void synchr_alpha_beta_to_theta_phi();
void synchr_theta_phi_to_alpha_beta();

public:

//qubit(); //constructeur par défaut
qubit (complexe alpha_ , complexe beta_);

complexe get_alpha() ;
complexe get_beta() ;
double get_abs_alpha() ;
double get_abs_beta() ;

void set_alpha(complexe alpha_);
void set_beta(complexe beta_);

void display();


 };

#endif