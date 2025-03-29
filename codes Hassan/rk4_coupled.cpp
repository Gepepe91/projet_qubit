#include"qubit.h"
#include"matrice.h"
#include<complex>
#include<cmath>
#include<vector>
#include<fstream>

using complexe = std::complex<double>;

//déclaration des fonctions
complexe f(const double& t , const complexe& alpha , const complexe& beta);
complexe g(const double& t , const complexe& alpha , const complexe& beta);
double omega_0 = 1.76e9;
double omega_1 = 1.e9;
double omega = 1.76e9;
int n = 1000;
double dt = 1.e-9;

void rk4_couple(std::vector<complexe>& alpha, std::vector<complexe>& beta , std::vector<double>& t , const double& dt , const int& n);



int main(){


    qubit q(2*M_PI/3,0); // qubit préparé dans (|0> + |1>)/sqrt(2)

    // dx/dt = f(x,y,t)
    // dy/dt = g(x,y,t)


    
    std::vector<double> t(n) ;
    std::vector<complexe> alpha(n) , beta(n);

    //Conditions initiales
    t[0] = 0; //T_initial, c'est 0
    alpha[0] = q.get_alpha();
    beta[0] = q.get_beta();


    rk4_couple(alpha , beta , t , dt , n);

    

    return 0;
}

//définitions des fonctions
complexe f(const double& t , const complexe& alpha , const complexe& beta){
    return complexe (0,-1) * ((omega_0/2) * alpha + (omega_1/2) * exp(complexe(0,- omega*t)) * beta);
}
complexe g(const double& t , const complexe& alpha , const complexe& beta){
    return complexe (0,-1) * ((omega_1/2) * exp(complexe(0,omega*t)) * alpha - (omega_0/2) * beta);
}

void rk4_couple(std::vector<complexe>& alpha, std::vector<complexe>& beta , std::vector<double>& t , const double& dt , const int& n){
    complexe F1 , F2 , F3 , F4 , G1 , G2 , G3 , G4 = 0.;

    std::ofstream fichier;
    fichier.open("dynamic_data.csv");
    fichier << "temps alpha beta norm_alpha norm_beta" << std::endl;

    for (int i = 0 ; i<n-1 ; i++){

        fichier << t[i] << " " << alpha[i] << " " << beta[i] << " " << norm(alpha[i]) << " " << norm(beta[i]) <<  std::endl;

        F1 = f(t[i] , alpha[i] , beta[i]);
        G1 = g(t[i] , alpha[i] , beta[i]);

        F2 = f(t[i] + (dt/2.) , alpha[i] + F1*(dt/2.) , beta[i] + G1 * (dt/2.));
        G2 = g(t[i] + (dt/2.) , alpha[i] + F1*(dt/2.) , beta[i] + G1 * (dt/2.));

        F3 = f(t[i] + (dt/2.) , alpha[i] + F2*(dt/2.) , beta[i] + G2 * (dt/2.));
        G3 = g(t[i] + (dt/2.) , alpha[i] + F2*(dt/2.) , beta[i] + G2 * (dt/2.));

        F4 = f(t[i] + dt , alpha[i] + F3 * dt , beta[i] + G3*dt);
        G4 = g(t[i] + dt , alpha[i] + F3 * dt , beta[i] + G3*dt);

        alpha[i+1] = alpha[i] + (dt/6.)*(F1 + 2.*F2 + 2.*F3 + F4);
        beta[i+1] = beta[i] + (dt/6.)*(G1 + 2.*G2 + 2.*G3 + G4);
        t[i+1] = t[i] + dt;

    }

    fichier.close();
}