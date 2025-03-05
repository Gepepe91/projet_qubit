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
    fichier.open("intrique.csv") ;

    fichier << "p N0 N1 N2 N3 estimation xi_exact ECM exact_ECM fisher fisher_exact" << std::endl;

    for (double p=1e-4 ; p < 0.99 ; p = p + 0.01) {
	    int N0 = 0 ;
        int N1 = 0 ;
        int N2 = 0 ;
        int N3 = 0 ;
        double estimation, ECM, ECM_exact, fisher, fisher_exact = 0. ;
	    double proba_0 = 1./2 * (1 - p/2 + (1-p) * std::cos (xi)) ;
        double proba_1 = p/4 ;
        double proba_2 = 1./2 * (1 - p/2 - (1-p) * std::cos (xi)) ;
        double proba_3 = p/4 ;
        double somme = proba_0 + proba_1 + proba_2 + proba_3 ;

	    for (int i=0 ; i < N ; i++) {
            double alea = dis(gen) ;
		    if (alea <= proba_0) {
			    N0 += 1 ;
            }
            else if ((alea <= proba_1 + proba_0) && (proba_0 <= alea)) {
                N1 += 1 ;
            }
            else if ((alea <= somme - proba_3) && (proba_1 + proba_0 <= alea)) {
                N2 += 1 ;
            }
            else if ((alea <= somme) && (somme - proba_3 <= alea)) {
                N3 += 1 ;
            }
        }

        double A = 1-p ;
        double B = 1 - p/2 ;
        estimation = acos((B*(N0-N2)) / (A*(N0 + N2))) ; 
        ECM = (B*B - A*A*cos(estimation)*cos(estimation)) / (B*A*A*sin(estimation)*sin(estimation)*N) ;
        ECM_exact = (B*B - A*A*cos(xi)*cos(xi)) / (B*A*A*sin(xi)*sin(xi)*N) ;
        fisher = (B*A*A*sin(estimation)*sin(estimation)) / (B*B - A*A*cos(estimation)*cos(estimation)) ;
        fisher_exact = (B*A*A*sin(xi)*sin(xi)) / (B*B - A*A*cos(xi)*cos(xi)) ;

        fichier << p << " " << N0 << " " << N1 << " " << N2 << " " << N3 << " " << estimation << " " << xi << " " << ECM << " " << ECM_exact << " " << fisher << " " << fisher_exact << std::endl ;
    }
    fichier.close() ;

    return 0 ;
    
}