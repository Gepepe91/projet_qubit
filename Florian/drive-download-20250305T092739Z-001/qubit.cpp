#include <iostream>
#include <random>
#include <fstream>

int main() {
    std::random_device rd;   
    std::mt19937 gen(rd());  // Générateur Mersenne Twister
    std::uniform_real_distribution<double> dis(0.0, 1.0); // Distribution uniforme entre 0 et 1

    int N = 200 ;
    double xi = M_PI/4 ;

    std::ofstream fichier ;
    fichier.open("qubit.csv") ;

    fichier << "p N0 estimation exact ECM exact_ECM fisher fisher_exact" << std::endl;
    
    for (double p=1e-4 ; p < 0.99 ; p = p + 0.01) {
	    int N0 = 0 ;
	    double projection = 1./2 * (1 + (1 - p) * std::cos (xi)) ;

	    for (int i=0 ; i < N ; i++) {
		    if (dis(gen) <= projection) {
			    N0 += 1 ;
            }
        }

        double A = 1-p ;
        double estimation, exact, ECM, exact_ECM, fisher, fisher_exact = 0.0 ;
	    estimation = std::acos((2*N0-N)/(N*A)) ;
        exact = xi ;
        ECM = (1 - A*A*cos(estimation)*cos(estimation))/(A*A*sin(estimation)*sin(estimation)*N) ;
        exact_ECM = (1 - A*A*cos(xi)*cos(xi))/(A*A*sin(xi)*sin(xi)*N) ;
        fisher = (A*A*sin(estimation)*sin(estimation))/ (1 - A*A*cos(estimation)*cos(estimation)) ;
        fisher_exact = (A*A*sin(xi)*sin(xi))/ (1 - A*A*cos(xi)*cos(xi)) ;
        fichier << p  << "," << N0 << "," << estimation << "," << exact << "," << ECM << "," << exact_ECM << "," << fisher << "," << fisher_exact << std::endl;
    }
    fichier.close() ;

    return 0 ;
    
}