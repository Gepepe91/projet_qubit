#ifndef QUBIT_H
#define QUBIT_H

#include <complex> //utilisation de la classe complex<double> de std
#include"matrice.h"

class qubit

/*
un qubit peut se modéliser de cette façon :
|psy(t)> = alpha|0> + beta|1> , alpha et beta std::complex<double>s --> utilisation de la classe std::complex<double>s créé en TD

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

std::complex<double> alpha , beta;
double theta , phi;

void synchr_alpha_beta_to_theta_phi();
void synchr_theta_phi_to_alpha_beta();

public:

qubit(); //constructeur par défaut
qubit (std::complex<double> alpha_ , std::complex<double> beta_);

std::complex<double> get_alpha() ;
std::complex<double> get_beta() ;
double get_abs_alpha() ;
double get_abs_beta() ;

void set_alpha(std::complex<double> alpha_);
void set_beta(std::complex<double> beta_);

void display();

void transform(matrice m);


 };

#endif