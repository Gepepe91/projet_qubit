#include<fstream>

// utilisation de l'équation de Schrödinger
// (d|Psi>)/(dt) = \frac{-i}{hbar} * H |Psi>
qubit compute_derivative(qubit q, matrice H) {
    // Appliquer la matrice H au qubit
    qubit new_q = H*q;
    // multiplication par -i / hbar
    new_q = new_q * (complexe (0,-1) * (1/hbar));
    return new_q;
}


//évolution par rk4 d'un qubit
 void rk4(qubit& q, matrice H, double dt) {
    //initialisation à theta = 0 et phi = 0
    qubit k1 , k2 , k3 , k4 ;
    // k1
    k1 = compute_derivative(q, H);
    // k2
    k2 = compute_derivative(q + k1*(0.5*dt), H);
    // k3
    k3 = compute_derivative(q + k2*(0.5*dt) , H);
    // k4
    k4 = compute_derivative(q + k3*dt , H);
    // Mise à jour de l'état du qubit
    q = q + (k1 + k2*2 + k3*2 + k4) * (dt/6); // l'ordre est un peu chiant, difficulté de la surcharge en C++
}


// simulation d'un qubit dans un champ statique
void simulation_static(qubit q , matrice H , double dt , double T_max){
    //Ouverture du fichier
    std::ofstream fichier;
    fichier.open("data.csv");

    //boucle for, simulation à proprement parlé
    fichier << "Time theta phi abs_alpha2 abs_beta2" << std::endl;
    for (double t = 0; t <= T_max; t += dt) {
        fichier << t  << " " << q.get_theta() << " " << q.get_phi() << " " << q.get_abs_alpha2() << " " << q.get_abs_beta2() << std::endl;
        rk4(q, H, dt); // Mise à jour de |psi> = notre qubit avec RK4
        q.normalize(); // renormalisation, sinon le qubit se dénormalise au cours du temps (erreur numérique)
    };
    //Fermeture du fichier
    fichier.close();
}









// Définition des équations différentielles couplées
// utilisées dans simulation_dynamic
complexe f(double t, complexe alpha, complexe beta, double omega , double omega_0 , double omega_1) {
    return complexe(0, -1) * ((omega_0 / 2) * alpha + (omega_1 / 2) * exp(complexe(0, -omega * t)) * beta);
}

complexe g(double t, complexe alpha, complexe beta, double omega , double omega_0 , double omega_1) {
    return complexe(0, -1) * ((omega_1 / 2) * exp(complexe(0, omega * t)) * alpha - (omega_0 / 2) * beta);
}






void simulation_dynamic(qubit q_init, double omega , double omega_0 , double omega_1, double dt, double T_init, double T_max, int n) {
    std::vector<double> t(n);
    std::vector<complexe> alpha(n), beta(n);

    // Conditions initiales
    t[0] = T_init ;
    alpha[0] = q_init.get_alpha();
    beta[0] = q_init.get_beta();



    std::ofstream fichier("dynamic_data.csv");
    fichier << "Time alpha beta abs_alpha2 abs_beta2" << std::endl;

    for (int i = 0; i < n - 1; i++) {
        fichier << t[i] << " " << alpha[i] << " " << beta[i] << " " 
                << norm(alpha[i]) << " " << norm(beta[i]) << std::endl;

        complexe F1 = f(t[i], alpha[i], beta[i], omega , omega_0 , omega_1);
        complexe G1 = g(t[i], alpha[i], beta[i], omega , omega_0 , omega_1);

        complexe F2 = f(t[i] + dt / 2.0, alpha[i] + F1 * (dt / 2.0), beta[i] + G1 * (dt / 2.0), omega , omega_0 , omega_1);
        complexe G2 = g(t[i] + dt / 2.0, alpha[i] + F1 * (dt / 2.0), beta[i] + G1 * (dt / 2.0),  omega , omega_0 , omega_1);

        complexe F3 = f(t[i] + dt / 2.0, alpha[i] + F2 * (dt / 2.0), beta[i] + G2 * (dt / 2.0), omega , omega_0 , omega_1);
        complexe G3 = g(t[i] + dt / 2.0, alpha[i] + F2 * (dt / 2.0), beta[i] + G2 * (dt / 2.0), omega , omega_0 , omega_1);

        complexe F4 = f(t[i] + dt, alpha[i] + F3 * dt, beta[i] + G3 * dt, omega , omega_0 , omega_1);
        complexe G4 = g(t[i] + dt, alpha[i] + F3 * dt, beta[i] + G3 * dt, omega , omega_0 , omega_1);

        alpha[i + 1] = alpha[i] + (dt / 6.0) * (F1 + 2.0 * F2 + 2.0 * F3 + F4);
        beta[i + 1] = beta[i] + (dt / 6.0) * (G1 + 2.0 * G2 + 2.0 * G3 + G4);
        t[i + 1] = t[i] + dt;
    }

    fichier.close();
}








// fonction de préparation d'un qubit dans l'état souhaité
// très similaire à la simulation en champ oscillant
// retourne un qubit
qubit preparation_etat_plus(qubit q_init, double omega , double omega_0 , double omega_1, double dt, double T_init, double T_max, int n) {
    std::vector<double> t(n);
    std::vector<complexe> alpha(n), beta(n);
    // On a du gâchis d'allocation de mémoire, la majorité des cases mémoires de nos vecteurs ne vont pas être
    // utilisées

    // Conditions initiales
    t[0] = T_init ;
    alpha[0] = q_init.get_alpha();
    beta[0] = q_init.get_beta();

    int rang=0; // pour obtenir les dernière valeurs de alpha et beta de l'itération



    std::ofstream fichier("preparation_qubit.csv");
    fichier << "Time alpha beta abs_alpha2 abs_beta2" << std::endl;

        
        for (int i = 0; i < n - 1; i++) {
            fichier << t[i] << " " << alpha[i] << " " << beta[i] << " " 
                    << norm(alpha[i]) << " " << norm(beta[i]) << std::endl;

            if (0.5 - 1e-6 <= norm(alpha[i]) && norm(alpha[i]) <= 0.5 + 1e-6){ //la précision ici est arbitraire
                rang = i;
                break;
            };

        /*
            if (0.5 - 1e-4 <= norm(alpha[i]) && norm(alpha[i]) <= 0.5 + 1e-4) {
                double epsilon = abs(norm(alpha[i]) - norm(alpha[i - 1]));
                if (0.5 - epsilon <= norm(alpha[i]) && norm(alpha[i]) <= 0.5 + epsilon) {
                    break;
                }
        */
                
                
                /*
                deux if pour break la boucle d'itération
                On exploite le fait que proche de norm(alpha)=norm(beta)=0.5 , nos courbes sont linéaires,
                donc la précision du pas est constante
                1. if P1 dans [0.5 - 1e-4 ; 0.5 + 1e4]
                        epsilon = norm(alpha[i] - norm(alpha[i-1]))
                        if P1 dans [0.5 - epsilon ; 0.5 + epsilon]
                            break !
                */
        
        

        complexe F1 = f(t[i], alpha[i], beta[i], omega , omega_0 , omega_1);
        complexe G1 = g(t[i], alpha[i], beta[i], omega , omega_0 , omega_1);

        complexe F2 = f(t[i] + dt / 2.0, alpha[i] + F1 * (dt / 2.0), beta[i] + G1 * (dt / 2.0), omega , omega_0 , omega_1);
        complexe G2 = g(t[i] + dt / 2.0, alpha[i] + F1 * (dt / 2.0), beta[i] + G1 * (dt / 2.0),  omega , omega_0 , omega_1);

        complexe F3 = f(t[i] + dt / 2.0, alpha[i] + F2 * (dt / 2.0), beta[i] + G2 * (dt / 2.0), omega , omega_0 , omega_1);
        complexe G3 = g(t[i] + dt / 2.0, alpha[i] + F2 * (dt / 2.0), beta[i] + G2 * (dt / 2.0), omega , omega_0 , omega_1);

        complexe F4 = f(t[i] + dt, alpha[i] + F3 * dt, beta[i] + G3 * dt, omega , omega_0 , omega_1);
        complexe G4 = g(t[i] + dt, alpha[i] + F3 * dt, beta[i] + G3 * dt, omega , omega_0 , omega_1);

        alpha[i + 1] = alpha[i] + (dt / 6.0) * (F1 + 2.0 * F2 + 2.0 * F3 + F4);
        beta[i + 1] = beta[i] + (dt / 6.0) * (G1 + 2.0 * G2 + 2.0 * G3 + G4);
        t[i + 1] = t[i] + dt;
    }

    fichier.close();

    q_init.set_alpha(alpha[rang]) ; q_init.set_beta(beta[rang]);
    return q_init;
}






// estimation Fisher avec un seul qubit
void one_qubit(int N , double xi) {

    //int N ; nb d'itération de la mesure
    //double xi ; valeur excate du déphasage

    //génération random
    std::random_device rd;   
    std::mt19937 gen(rd());  // Générateur Mersenne Twister
    std::uniform_real_distribution<double> dis(0.0, 1.0); // Distribution uniforme entre 0 et 1

    //ouverture fichier et écriture de la première ligne
    std::ofstream fichier ;
    fichier.open("estimation.csv") ;
    fichier << "p N0 estimation exact ECM exact_ECM fisher fisher_exact" << std::endl;

    // boucle sur différentes probabilités
    // l'utilisation de la fonction arccos oblige d'aller de 0 + epsilon à 1 - epsilon
    // p caractérise à quel point notre bruit peut affecter le qubit déphasé
    for (double p=1e-4 ; p < 0.99 ; p = p + 0.01) {

        // initialisation du compteur du nb de fois où on projette sur |0>
	    int N0 = 0 ;

        //proba de mesurer la valeur propre associé à "|+>" , l'état préparé
	    double proba_0 = 1./2 * (1 + (1 - p) * std::cos (xi)) ;

        // ajout du nombre de fois où on mesure valeur propre de |+>
        // si notre  nb random est inférieur à la proba, on considère que nous mesurons |+> ==> N0 += 1;
	    for (int i=0 ; i < N ; i++) {
		    if (dis(gen) <= proba_0) {
			    N0 += 1 ;
            }
        }

        // A est simplement une variable intermédiaire pour améliorer la lisibilité des calculs
        double A = 1-p ;
        //iniitialisation des varaibles
        double estimation_dephasage, exact, ECM, exact_ECM, fisher, fisher_exact = 0.0 ;

        //calculs, issus de l'article :

        //estimation du paramètre de déphasage
	    estimation_dephasage = std::acos((2*N0-N)/(N*A)) ;
        // déphasage exact
        exact = xi ;
        //Ecart quadratique moyen avec l'estimation du déphasage
        ECM = (1 - A*A*cos(estimation_dephasage)*cos(estimation_dephasage))/(A*A*sin(estimation_dephasage)*sin(estimation_dephasage)*N) ;
        //ecart quadratique moyen avec la valeur exacte du déphasage
        exact_ECM = (1 - A*A*cos(xi)*cos(xi))/(A*A*sin(xi)*sin(xi)*N) ;
        //information de Fisher pour le paramètre de déphasage estimé
        fisher = (A*A*sin(estimation_dephasage)*sin(estimation_dephasage))/ (1 - A*A*cos(estimation_dephasage)*cos(estimation_dephasage)) ;
        //information de Fisher pour le paramètre de déphasage exact
        fisher_exact = (A*A*sin(xi)*sin(xi))/ (1 - A*A*cos(xi)*cos(xi)) ;

        if (!std::isnan(estimation_dephasage)){
            //écriture dans le fichier en excluant les nan
            fichier << p  << "," << N0 << "," << estimation_dephasage << "," << exact << "," << ECM << "," << exact_ECM << "," << fisher << "," << fisher_exact << std::endl;
        }
    }


    fichier.close() ;
    
}

// estimation Fisher avec 2 qubits
void two_qubits(int N , double xi) {

    //int N ; nb d'itération de la mesure
    //double xi ; valeur excate du déphasage

    // random entre 0 et 1
    std::random_device rd;   
    std::mt19937 gen(rd());  // Générateur Mersenne Twister
    std::uniform_real_distribution<double> dis(0.0, 1.0); // Distribution uniforme entre 0 et 1

    std::ofstream fichier ;
    fichier.open("intrique.csv") ;

    // 2 qubit intriqués --> espace de Hilbert de dimension 2^2=4
    // Plus possible de discriminer simplement

    fichier << "p N0 N1 N2 N3 estimation xi_exact ECM exact_ECM fisher fisher_exact" << std::endl;

    for (double p=1e-4 ; p < 0.99 ; p = p + 0.01) { //boucle de p presque 0 à 1, arccos...
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

        //initialisation
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

        // A et B pour améliorer la lisibilité des calculs
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

        fichier << p << " " << N0 << " " << N1 << " " << N2 << " " << N3 << " " << estimation << " " << xi << " " << ECM << " " << ECM_exact << " " << fisher << " " << fisher_exact << std::endl ;
    }
    fichier.close() ;
    
}