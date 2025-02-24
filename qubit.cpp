#include<iostream>
#include<cmath>
#include<complex>
#include"qubit.h"

using complex = std::complex<double>;

// ATTENTION : abs et arg sont des fonctions qui s'appliquent aux complexes, pas des méthodes


qubit::qubit() : alpha(complex(0,0)) , beta(complex(0,0)){} //constructeur par défaut
qubit::qubit(complex alpha_ , complex beta_) : alpha(alpha_) , beta(beta_){}

void qubit::synchr_alpha_beta_to_theta_phi(){
    theta = 2*acos(abs(alpha));
    phi = arg(beta) - arg(alpha);
}
void qubit::synchr_theta_phi_to_alpha_beta(){
    alpha = cos(theta/2);
    beta = exp(complex(0,phi)) * sin(theta/2);
}

void qubit::set_alpha(complex alpha_){
    alpha = alpha_;
    synchr_alpha_beta_to_theta_phi();
}

void qubit::set_beta(complex beta_){
    beta = beta_;
    synchr_alpha_beta_to_theta_phi();
}

complex qubit::get_alpha() {
    return alpha;
}

complex qubit::get_beta() {
    return beta;
}

double qubit::get_theta() {
    synchr_alpha_beta_to_theta_phi();
    return theta;
}
double qubit::get_phi() {
    synchr_alpha_beta_to_theta_phi();
    return phi;
}

double qubit::get_abs_alpha() {
    return abs(alpha);
}

double qubit::get_abs_beta() {
    return abs(beta);
}


void qubit::display(){
    std::cout << alpha; std::cout << beta;
}

void qubit::transform(matrice m){
    alpha = m.get_element(0,0)*alpha + m.get_element(0,1)*beta;
    beta = m.get_element(1,0)*alpha + m.get_element(1,1)*beta;
}

void qubit::normalize() {
    double constante_norm = std::sqrt(std::norm(alpha) + std::norm(beta));
    alpha /= constante_norm;
    beta /= constante_norm;
}

// Surcharge de l'opérateur * pour multiplier par un scalaire complexe
qubit qubit::operator*(const std::complex<double>& scalaire) {
    // Multiplie chaque composant par le scalaire
    return qubit(alpha * scalaire, beta * scalaire);
}
// Surcharge de l'opérateur * pour multiplier par un scalaire complexe
qubit qubit::operator*(const double& scalaire){
    return qubit(alpha * scalaire , beta * scalaire);
}
//Surcharge de l'opérateur + entre qubit
qubit qubit:: operator+(const qubit& r){
    return qubit(alpha + r.alpha , beta + r.beta);
}