#include <iostream>
#include <random>
#include <fstream>



// POur 1 qubit

int main() {

    //génération random
    std::random_device rd;   
    std::mt19937 gen(rd());  // Générateur Mersenne Twister
    std::uniform_real_distribution<double> dis(0.0, 1.0); // Distribution uniforme entre 0 et 1

    int N = 200 ; // nb d'itération de la mesure
    double xi = M_PI/4 ; //valeur excate du déphasage

    std::ofstream fichier ;
    fichier.open("estimation.csv") ;

    fichier << "p N0 estimation exact ECM exact_ECM fisher fisher_exact" << std::endl;


    //variation proba de presque 0 à presque 1 pour les fonctions arccos
    // p caractérise à quel point notre bruit peut affcter le qubit déphasé
    for (double p=1e-4 ; p < 0.99 ; p = p + 0.01) {
	    int N0 = 0 ; //compteur nb de projection sur état |0>

        //proba de mesurer vp ass. à psi_0 l'état préparé, ici "|+>"
	    double proba_0 = 1./2 * (1 + (1 - p) * std::cos (xi)) ;

	    for (int i=0 ; i < N ; i++) {
		    if (dis(gen) <= proba_0) { //si notre nb random est inf à la proba --> mesure vp de |+> et on ajoute +1 à N0
			    N0 += 1 ;
            }
        }

        double A = 1-p ;
        double estimation, exact, ECM, exact_ECM, fisher, fisher_exact = 0.0 ;
	    estimation = std::acos((2*N0-N)/(N*A)) ; //estimation paramètre de déphasage
        exact = xi ; //donnée initialement, connue
        //Ecart quadratique moyen avec l'estimation du déphasage
        ECM = (1 - A*A*cos(estimation)*cos(estimation))/(A*A*sin(estimation)*sin(estimation)*N) ;
        //ecart quadratique moyen avec la valeur exacte du déphasage
        exact_ECM = (1 - A*A*cos(xi)*cos(xi))/(A*A*sin(xi)*sin(xi)*N) ;
        //info de fisher pour le param estimé
        fisher = (A*A*sin(estimation)*sin(estimation))/ (1 - A*A*cos(estimation)*cos(estimation)) ;
        // me^me chose mais param exacte
        fisher_exact = (A*A*sin(xi)*sin(xi))/ (1 - A*A*cos(xi)*cos(xi)) ;
        fichier << p  << "," << N0 << "," << estimation << "," << exact << "," << ECM << "," << exact_ECM << "," << fisher << "," << fisher_exact << std::endl;
    }
    fichier.close() ;

    return 0 ;
    
}