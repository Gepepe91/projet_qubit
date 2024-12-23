#include<iostream>
#include<cmath>
#include"qubit.h"
#include"complexe.h"
#include"surcharge_complexe.h"

qubit::qubit(complexe alpha_ , complexe beta_) : alpha(alpha_) , beta(beta_){}

void qubit::synchr_alpha_beta_to_theta_phi(){
    theta = 2*acos(alpha.get_norm());
    phi = beta.get_phase() - alpha.get_phase();
}
void qubit::synchr_theta_phi_to_alpha_beta(){
}

void qubit::set_alpha(complexe alpha_){
    alpha = alpha_;
    synchr_alpha_beta_to_theta_phi();
}

void qubit::set_beta(complexe beta_){
    beta = beta_;
    synchr_alpha_beta_to_theta_phi();
}

complexe qubit::get_alpha() {
    return alpha;
}

complexe qubit::get_beta() {
    return beta;
}

double qubit::get_abs_alpha() {
    return alpha.get_part_reelle();
}

double qubit::get_abs_beta() {
    return beta.get_part_imaginaire();
}

void qubit::display(){
    alpha.afficher(); beta.afficher();
}

void qubit::transform(matrice m){
    alpha = m.get_element(0,0)*alpha + m.get_element(0,1)*beta;
    beta = m.get_element(1,0)*alpha + m.get_element(1,1)*beta;
}