#include <iostream>
#include <random>
#include <fstream>
#include <cmath>
#include <chrono>

// Pour 1 qubit

int main() {
    
    // random entre 0 et 1
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); //génération random en fonction de l'heure
    std::mt19937 gen(seed); //Générateur Mersenne Twister 
    std::uniform_real_distribution<double> dis(0.0, 1.0); // Distribution uniforme entre 0 et 1

    int N; // nb d'itération de la mesure
    std::cout << "Entrez le nombre de qubits préparés et mesurés N :"; std::cin >> N;
    double xi = M_PI/4 ; //valeur exacte du déphasage

    std::ofstream fichier ;
    fichier.open("estimation_single.csv") ; // fichier qui récupère les données pour 

    fichier << "p N0 estimation exact ECM exact_ECM fisher fisher_exact" << std::endl;


    //variation proba. de 0 à 1 EXCLUS pour éviter les divergences avec arccos
    // p caractérise à quel point notre bruit peut affcter le qubit déphasé
    for (double p=1e-4 ; p < 0.99 ; p = p + 0.005) {
	    int N0 = 0 ; //compteur nb de projection sur état |0>

        //proba de mesurer la valeur propre associée à psi_0 l'état préparé, ici "|+>"
	    double proba_0 = 1./2 * (1 + (1 - p) * std::cos(xi)) ;

	    for (int i=0 ; i < N ; i++) {
		    if (dis(gen) <= proba_0) { //si notre nb random est inférieur à la proba --> mesure valeur propre de |+> et on ajoute +1 à N0
			    N0 += 1 ;
            }
        }

        double A = 1-p ;
        double estimation, exact, ECM, exact_ECM, fisher, fisher_exact = 0.0 ;
	    estimation = std::acos((2*N0-N)/(N*A)) ; //estimation paramètre de déphasage
        exact = xi ; //déphasage xi exact donnée initialement, connue

        //Ecart quadratique moyen avec l'estimation du déphasage xi
        ECM = (1 - A*A*cos(estimation)*cos(estimation))/(A*A*sin(estimation)*sin(estimation)*N) ;

        //ecart quadratique moyen avec la valeur exacte du déphasage
        exact_ECM = (1 - A*A*cos(xi)*cos(xi))/(A*A*sin(xi)*sin(xi)*N) ;

        //info de fisher avec l'estimation du déphasage xi
        fisher = (A*A*sin(estimation)*sin(estimation))/ (1 - A*A*cos(estimation)*cos(estimation)) ;

        //info de fisher avec la valeur exacte du déphasage xi 
        fisher_exact = (A*A*sin(xi)*sin(xi))/ (1 - A*A*cos(xi)* cos(xi)) ;

        // Condition pour ne pas considérer les valeurs NaN 
        if (!std::isnan(estimation)) {
            fichier << p  << "," << N0 << "," << estimation << "," << exact << "," << ECM << "," << exact_ECM << "," << fisher << "," << fisher_exact << std::endl;
        }
    }
    fichier.close();

    return 0 ;
}
