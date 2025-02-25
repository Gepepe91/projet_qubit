#ifndef qubit_H
#define qubit_H

#include <complex>

class qubit {
private:

//Attributs caractérisant notre qubit + fonctions de synchronisation

    std::complex<double> alpha , beta;  // Amplitudes pour |0> et |1>

    double theta , phi;          // Paramètres theta et phi

    void synchr_alpha_beta_to_theta_phi();
    void synchr_theta_phi_to_alpha_beta();

public:

    // Constructeurs

    // Constructeur par défaut
    qubit();
    // Constructeur pour initialiser un qubit avec des paramètres θ et φ
    qubit(double theta, double phi);
    // Constructeur pour initialiser avec alpha et beta
    qubit(std::complex<double> alpha , std::complex<double> beta);

    //Getters

    // Getter pour obtenir l'amplitude de |0> et |1>
    complex<double> get_alpha();
    complex<double> get_beta();
    //Getter pour obtenir les angles theta et phi
    double get_theta();
    double get_phi();
    //Getter pour obtenir les probabilités de |0> et |1>
    double get_abs_alpha2();
    double get_abs_beta2();

    //Setter

    //setter de alpha et beta
    void set_alpha(std::complex<double> alpha_);
    void set_beta(std::complex<double> beta_);
    //setter de theta et phi
    void set_theta(double theta_);
    void set_phi(double phi_);

    //Méthodes de calculs et modifications
    void tranform(matrice m); //matrice étant une classe matrice faite maison
    void normalize(); //Permet d'éviter l'accumulation d'erreurs numériques sur notre qubit et s'assurer qu'il reste normé

    //Méthodes d'affichages
    void display();

    //Surcharges d'opérateurs
    qubit operator*(const std::complex<double>& scalaire);
    qubit operator*(const double& scalaire);
    qubit operator+(const qubit& r);

};

#endif