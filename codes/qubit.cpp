#include<iostream>
#include<cmath>
#include "qubit.h"
//#include"matrice.h"

using complex = std::complex<double>;

//Constructeurs

//par défaut
qubit::qubit() : theta(0) , phi(0){} //Par défaut, les angles sont nuls, pôle nord dans le plan 0x
//Par alpha et beta
qubit::qubit(complex alpha_ , complex beta_) : alpha(alpha_) , beta(beta_){
    synchr_alpha_beta_to_theta_phi();
}
//Par theta et phi
qubit::qubit(double theta_ , double phi_) : theta(theta_) , phi(phi_){
    synchr_theta_phi_to_alpha_beta();
}

//Méthodes de synchronisation
void qubit::synchr_alpha_beta_to_theta_phi(){
    if (abs(alpha) >= 1 - 1e-12) {
        theta = 0.;
        phi = 0;}
    else if (abs(beta) >= 1 - 1e-12){
        theta = M_PI;
        phi = 0;
    }
    else{
        theta = 2*acos(abs(alpha));}
    phi = arg(beta) - arg(alpha);
}
void qubit::synchr_theta_phi_to_alpha_beta(){
    alpha = cos(theta/2);
    beta = exp(complex(0,phi)) * sin(theta/2);
}

//Getters

//alpha et beta
complex qubit::get_alpha() {
    synchr_theta_phi_to_alpha_beta();
    return alpha;
}
complex qubit::get_beta() {
    synchr_theta_phi_to_alpha_beta();
    return beta;
}
//theta et phi
double qubit::get_theta() {
    synchr_alpha_beta_to_theta_phi();
    return theta;
}
double qubit::get_phi() {
    synchr_alpha_beta_to_theta_phi();
    return phi;
}
//probabilités de |0> et |1>
double qubit::get_abs_alpha2() {
    synchr_theta_phi_to_alpha_beta();
    return norm(alpha);
}
double qubit::get_abs_beta2() {
    synchr_theta_phi_to_alpha_beta();
    return norm(beta);
}

//Setters

//alpha et beta
void qubit::set_alpha(complex alpha_){
    alpha = alpha_;
    synchr_alpha_beta_to_theta_phi();
}
void qubit::set_beta(complex beta_){
    beta = beta_;
    synchr_alpha_beta_to_theta_phi();
}
//theta et phi
void qubit::set_theta(double theta_){
    theta = theta_;
    synchr_theta_phi_to_alpha_beta();
}
void qubit::set_phi(double phi_){
    phi = phi_;
    synchr_theta_phi_to_alpha_beta();
}

//méthodes de calculs et modifications


//Normalisation du qubit
void qubit::normalize() {
    double constante_norm = std::sqrt(std::norm(alpha) + std::norm(beta));
    alpha /= constante_norm;
    beta /= constante_norm;
}

//Méthodes visuelles
void qubit::display(){
    synchr_theta_phi_to_alpha_beta();
    std::cout << std::endl;
    std::cout << alpha << " coefficient complexe de |0>" << std::endl ;
    std::cout << beta << " coefficient complexe de |1>" << std::endl;
    std::cout << std::endl;
}
void qubit::display_angles(){
    synchr_alpha_beta_to_theta_phi();
    std::cout << std::endl;
    std::cout << theta << " angle theta" << std::endl;
    std::cout << phi << " angle phi" << std::endl;
    std::cout << std::endl;
}

//Surcharges d'opérateurs

// Surcharge de l'opérateur * pour multiplier par un scalaire complexe
qubit qubit::operator*(const complex& scalaire) {
    // Multiplie chaque composant par le scalaire
    return qubit(alpha * scalaire, beta * scalaire);
}
// Surcharge de l'opérateur * pour multiplier par un scalaire complexe
qubit qubit::operator*(const double& scalaire){
    return qubit(alpha * scalaire , beta * scalaire);
}
//Surcharge de l'opérateur + entre qubit
qubit qubit::operator+(const qubit& r){
    return qubit(alpha + r.alpha , beta + r.beta);
}

complex qubit::operator|(const qubit& q2){
    return conj(alpha) * q2.alpha + conj(beta) * q2.beta;
}
matrice qubit::operator&(const qubit& q2){
    return matrice (alpha * conj(q2.alpha) , alpha * conj(q2.beta) , beta * conj(q2.alpha) , beta * conj(q2.beta));
}

//matrice qubit::operator^(const qubit& q2){
//    return matrice (alpha * q2.alpha , alpha*q2.beta , beta*q2.alpha , beta*q2.beta);
//}
// ! A redéfinir, il faut obtenir une matrice 4*1