#include<iostream>
#include<cmath>
#include<complex>
#include"qubit.h"

using complex = std::complex<double>;

// ATTENTION : abs et arg sont des fonctions qui s'appliquent aux complexes, pas des méthodes


qubit::qubit(complex alpha_ , complex beta_) : alpha(alpha_) , beta(beta_){}

void qubit::synchr_alpha_beta_to_theta_phi(){
    theta = 2*acos(abs(alpha));
    phi = arg(beta) - arg(alpha);
}
void qubit::synchr_theta_phi_to_alpha_beta(){
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

double qubit::get_abs_alpha() {
    return alpha.real();
}

double qubit::get_abs_beta() {
    return beta.imag();
}

void qubit::display(){
    std::cout << alpha; std::cout << beta;
}

void qubit::transform(matrice m){
    alpha = m.get_element(0,0)*alpha + m.get_element(0,1)*beta;
    beta = m.get_element(1,0)*alpha + m.get_element(1,1)*beta;
}