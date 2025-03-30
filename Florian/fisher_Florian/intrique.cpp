#include <iostream>
#include <random>
#include <fstream>
#include <chrono>

// 2 qubit intriqués

int main() {
    // random entre 0 et 1
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); //génération random en fonction de l'heure
    std::mt19937 gen(seed); //Générateur Mersenne Twister 
    std::uniform_real_distribution<double> dis(0.0, 1.0); // Distribution uniforme entre 0 et 1

    int N = 200 ;
    double xi = M_PI/4 ;//valeur exacte déphasage

    std::ofstream fichier ;
    fichier.open("intrique.csv") ;

    // 2 qubit intriqués --> espace de Hilbert de dimension 2^2=4
    // Plus possible de discriminer simplement

    fichier << "p N0 N1 N2 N3 estimation xi_exact ECM exact_ECM fisher fisher_exact" << std::endl;

    for (double p=1e-4 ; p < 0.99 ; p = p + 0.005) { //boucle de p presque 0 à 1, arccos...
        //4 compteur pour les 4 états sur lesquels on projette
	    int N0 = 0 ;
        int N1 = 0 ;
        int N2 = 0 ;
        int N3 = 0 ;
        //proba pour les 4 mesures, |++> , |+-> , |-+> , |-->
	    double proba_0 = 1./2 * (1 - p/2 + (1-p) * std::cos (xi)) ;
        double proba_1 = p/4 ;
        double proba_2 = 1./2 * (1 - p/2 - (1-p) * std::cos (xi)) ;
        double proba_3 = p/4 ;
        // somme des probas, doit être 1
        double somme = proba_0 + proba_1 + proba_2 + proba_3 ;

        double estimation, ECM, ECM_exact, fisher, fisher_exact = 0. ;


        //problème, les probas des états se superposent, changement d'origine des probas, pas de discrim claire

        // 1 : mesure de P0
        // décalage de l'origine de P1 , 0 à P1 --> P0 à P0 + P1
        // etc avec P2 , P3 ...
        // Cela permet de discriminer les projections
        //normalement reste borné par 1

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

        estimation = acos((B*(N0-N2)) / (A*(N0 + N2))) ; //estimation du param xi
        //ecart quadra moyen du param estimé
        ECM = (B*B - A*A*cos(estimation)*cos(estimation)) / (B*A*A*sin(estimation)*sin(estimation)*N) ;
        //ecart quadra moyen du param exact
        ECM_exact = (B*B - A*A*cos(xi)*cos(xi)) / (B*A*A*sin(xi)*sin(xi)*N) ;
        // information de Fisher du param estimé
        fisher = (B*A*A*sin(estimation)*sin(estimation)) / (B*B - A*A*cos(estimation)*cos(estimation)) ;
        // info de Fisher du param exact
        fisher_exact = (B*A*A*sin(xi)*sin(xi)) / (B*B - A*A*cos(xi)*cos(xi)) ;

        if (!std::isnan(estimation)) {
            fichier << p << " " << N0 << " " << N1 << " " << N2 << " " << N3 << " " << estimation << " " << xi << " " << ECM << " " << ECM_exact << " " << fisher << " " << fisher_exact << std::endl ;
        }
    }
    fichier.close() ;

    return 0 ;
    
}