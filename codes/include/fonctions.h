#include<fstream>

/*
Nous avons créer ce fichier fonctions.h afin de rendre la navigation dans nos codes plus simple
En faisant #include"fonctions.h" au début de notre fichier main.cpp, c'est comme si nous concaténions
toutes les fonctions que nous utilisons ≠ des headers de nos classes.
*/

/*
1
Evolution d'un qubit seul non bruité dans un champ statique
*/

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





/*
2
Evolution d'un qubit seul non bruité dans un champ dynamique
*/



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







/*
Estimation Fisher
*/


// estimation Fisher avec un seul qubit
void one_qubit(int N , double xi) {

    //int N ; nb d'itération de la mesure
    //double xi ; valeur exacte du déphasage

    //génération random
    std::random_device rd;   
    std::mt19937 gen(rd());  // Générateur Mersenne Twister
    std::uniform_real_distribution<double> dis(0.0, 1.0); // Distribution uniforme entre 0 et 1

    //ouverture fichier et écriture de la première ligne
    std::ofstream fichier ;
    fichier.open("estimation.csv") ;
    fichier << "p N0 estimation exact ECM exact_ECM fisher fisher_exact" << std::endl;

    // boucle sur les différentes probabilités
    // l'utilisation de la fonction arccos oblige d'aller de 0 + epsilon à 1 - epsilon
    // p caractérise à quel point notre bruit peut affecter le qubit déphasé
    for (double p=1e-4 ; p < 0.99 ; p = p + 0.01) {

        // initialisation du compteur du nb de fois où on projette sur |0>
	    int N0 = 0 ;

        //proba de mesurer la valeur propre associé à "|+>" , l'état préparé
	    double proba_0 = 1./2 * (1 + (1 - p) * std::cos (xi)) ;

        // ajout du nombre de fois où on mesure la valeur propre de |+>
        // si notre  nb random est inférieur à la proba, on considère que nous mesurons |+> ==> N0 += 1;
	    for (int i=0 ; i < N ; i++) {
		    if (dis(gen) <= proba_0) {
			    N0 += 1 ;
            }
        }

        // A est une variable intermédiaire pour améliorer la lisibilité des calculs
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




// include des librairies d'aléatoire
static std::random_device rd;
static std::mt19937 gen(rd());


/*
Bruits 
*/

//Permet d'appeler les étapes RK4 couplées --> comme la simulation avec le champ dynamique
void rk4_step(complexe& alpha, complexe& beta, double t, double dt, double omega, double omega_0, double omega_1) {
    complexe F1 = f(t, alpha, beta, omega, omega_0, omega_1);
    complexe G1 = g(t, alpha, beta, omega, omega_0, omega_1);

    complexe F2 = f(t + dt / 2.0, alpha + F1 * (dt / 2.0), beta + G1 * (dt / 2.0), omega, omega_0, omega_1);
    complexe G2 = g(t + dt / 2.0, alpha + F1 * (dt / 2.0), beta + G1 * (dt / 2.0), omega, omega_0, omega_1);

    complexe F3 = f(t + dt / 2.0, alpha + F2 * (dt / 2.0), beta + G2 * (dt / 2.0), omega, omega_0, omega_1);
    complexe G3 = g(t + dt / 2.0, alpha + F2 * (dt / 2.0), beta + G2 * (dt / 2.0), omega, omega_0, omega_1);

    complexe F4 = f(t + dt, alpha + F3 * dt, beta + G3 * dt, omega, omega_0, omega_1);
    complexe G4 = g(t + dt, alpha + F3 * dt, beta + G3 * dt, omega, omega_0, omega_1);

    alpha = alpha + (dt / 6.0) * (F1 + 2.0 * F2 + 2.0 * F3 + F4);
    beta = beta + (dt / 6.0) * (G1 + 2.0 * G2 + 2.0 * G3 + G4);
}

// Fonction pour appliquer un bruit aléatoire au qubit
void appliquer_bruit_aleatoire(qubit& q, double niveau_bruit, double* bruittheta, double* bruitphi) {
    std::normal_distribution<double> dis(0, niveau_bruit);
    double bruit_theta = dis(gen)*M_PI/2;
    double bruit_phi = dis(gen)*M_PI;

    *bruittheta = bruit_theta;
    *bruitphi = bruit_phi;

    q.set_theta(q.get_theta() + bruit_theta);
    q.set_phi(q.get_phi() + bruit_phi);
    q.normalize();
}

void appliquer_bruit_aleatoire_ab(qubit& q, double niveau_bruit, double* bruitalpha, double* bruitbeta) {
    std::normal_distribution<double> dis(0, niveau_bruit);
    double bruit_alpha = dis(gen);
    double bruit_beta = dis(gen); 

    *bruitalpha = bruit_alpha;
    *bruitbeta = bruit_beta;

    // Appliquer un bruit sous forme d'une petite perturbation multiplicative
    q.set_alpha(q.get_alpha() * (1.0 + bruit_alpha));
    q.set_beta(q.get_beta() * (1.0 + bruit_beta));
    q.normalize();
}

// Fonction de simulation avec bruit

void simulation_avec_bruit(qubit q, matrice H, double dt, double T_init, double T_max, int n, double noise_level) {
    // Ouverture du fichier
    std::ofstream fichier;
    fichier.open("data_with_noise.csv");
    double t = T_init;
    double bruittheta = 0.0, bruitphi = 0.0;
    // Vérifiez si le fichier est ouvert correctement
    if (!fichier.is_open()) {
        std::cerr << "Erreur: Impossible d'ouvrir le fichier data_with_noise.csv" << std::endl;
        return;
    }

    // Boucle for, simulation à proprement parlé
    fichier << "Time theta phi abs_alpha2 abs_beta2" << std::endl;
    for (int i = 0; i < n - 1; i++) {
        fichier << t << " " << q.get_theta() << " " << q.get_phi() << " " << q.get_abs_alpha2() << " " << q.get_abs_beta2() << std::endl;
        //std::cout << "Écriture des données à t = " << t << ", theta = " << q.get_theta() << ", phi = " << q.get_phi() << ", abs_alpha2 = " << q.get_abs_alpha2() << ", abs_beta2 = " << q.get_abs_beta2() << std::endl; // Message de débogage
        rk4(q, H, dt); // Mise à jour de |psi> = notre qubit avec RK4
        appliquer_bruit_aleatoire(q, noise_level, &bruittheta, &bruitphi); // Appliquer le bruit aléatoire
        q.normalize(); // Renormalisation, sinon le qubit se dénormalise au cours du temps (erreur numérique)
        t = t + dt;
    }
    // Fermeture du fichier
    fichier.close();
}

// Fonction de simulation avec champ magnétique oscillant et bruit
void simulation_oscillating_magnetic_field_with_noise_and_correction( qubit q_init, double omega, double omega_0, 
    double omega_1, double dt, double T_init, double T_max, int n, double noise_level) {

    // Qubit bruité (expérience principale)
    qubit qubit_bruit = q_init;
    // Qubit corrigé
    qubit qubit_correction = q_init;

    double t = T_init;
    double bruit_theta = 0.0;
    double bruit_phi = 0.0;
    double ancien_bruit_theta = 0.0, ancien_bruit_phi = 0.0; // Bruit de l'itération précédente

    std::ofstream fichier("oscillating_magnetic_field_with_noise_and_correction.csv");

    // Check if the file is open
    if (!fichier.is_open()) {
        std::cerr << "Erreur: Impossible d'ouvrir le fichier oscillating_magnetic_field_with_noise.csv" << std::endl;
        return;
    }

    // Write header
    fichier << "Time abs_alpha_bruit abs_beta_bruit abs_alpha_corrige abs_beta_corrige bruit_alpha bruit_beta" << std::endl;

    // Simulation loop
    for (int i = 0; i < n - 1; i++) {
        // Write data to file
        fichier << t << " " << norm(qubit_bruit.get_alpha()) << " " << norm(qubit_bruit.get_beta()) << " "
                << norm(qubit_correction.get_alpha()) << " " << norm(qubit_correction.get_beta()) << " "
                << bruit_theta << " " << bruit_phi << std::endl;

        // Apply RK4 integration to evolve the state
        complexe alpha = qubit_bruit.get_alpha();
        complexe beta = qubit_bruit.get_beta();
        
        rk4_step(alpha,beta, t, dt, omega, omega_0, omega_1);

        // Update the noisy qubit
        qubit_bruit.set_alpha(alpha);
        qubit_bruit.set_beta(beta);
        qubit_bruit.normalize();

        // Apply noise to the qubit
        appliquer_bruit_aleatoire(qubit_bruit, noise_level, &bruit_theta, &bruit_phi);

        // Imprimer les valeurs de bruit pour le débogage
        std::cout << "Itération " << i << " | Bruit Alpha: " << bruit_theta << " | Bruit Beta: " << bruit_phi << std::endl;

        // Correct the noise
        double theta_corrige = qubit_correction.get_theta() - ancien_bruit_theta + bruit_theta;
        double phi_corrige = qubit_correction.get_phi() - ancien_bruit_phi+ bruit_phi;

        // Update the corrected qubit
        qubit_correction.set_theta(theta_corrige);
        qubit_correction.set_phi(phi_corrige);
        qubit_correction.normalize();

        complexe alpha_corr = qubit_correction.get_alpha();
        complexe beta_corr = qubit_correction.get_beta();

        rk4_step(alpha_corr, beta_corr , t, dt, omega, omega_0, omega_1);

        qubit_correction.set_alpha(alpha_corr);
        qubit_correction.set_beta(beta_corr);
        qubit_correction.normalize();

        // Mettre à jour les anciens bruits pour la prochaine itération
        ancien_bruit_theta = bruit_theta;
        ancien_bruit_phi = bruit_phi;

        // Update time
        t += dt;
    }

    // Close the file
    fichier.close();
}

void simulation_multiple_noise_levels(qubit q_init, double omega, double omega_0, double omega_1,
    double dt, double T_init, double T_max, int n) {
std::vector<double> noise_levels = {0.01, 0.02, 0.03, 0.04}; // Niveaux de bruit à tester
std::ofstream fichier("simulation_results.csv");

if (!fichier.is_open()) {
    std::cerr << "Erreur: Impossible d'ouvrir le fichier simulation_results.csv" << std::endl;
    return;
}

// Écrire l'en-tête du fichier CSV
fichier << "Noise Level,Time,abs_alpha_bruit,abs_beta_bruit,abs_alpha_corrige,abs_beta_corrige" << std::endl;

for (double noise_level : noise_levels) {
    qubit qubit_bruit = q_init;
    qubit qubit_correction = q_init;

    double t = T_init;
    double bruit_theta = 0.0;
    double bruit_phi = 0.0;
    double ancien_bruit_theta = 0.0, ancien_bruit_phi = 0.0;

    for (int i = 0; i < n - 1; i++) {
        fichier << noise_level << "," << t << ","
        << norm(qubit_bruit.get_alpha()) << "," << norm(qubit_bruit.get_beta()) << ","
        << norm(qubit_correction.get_alpha()) << "," << norm(qubit_correction.get_beta()) << std::endl;

        complexe alpha = qubit_bruit.get_alpha();
        complexe beta = qubit_bruit.get_beta();

        // Appliquer RK4 pour faire évoluer le qubit bruité
        rk4_step(alpha, beta, t, dt, omega, omega_0, omega_1);


        qubit_bruit.set_alpha(alpha);
        qubit_bruit.set_beta(beta);
        qubit_bruit.normalize();

        appliquer_bruit_aleatoire(qubit_bruit, noise_level, &bruit_theta, &bruit_phi);

        double theta_corrige = qubit_correction.get_theta() - ancien_bruit_theta + bruit_theta;
        double phi_corrige = qubit_correction.get_phi() - ancien_bruit_phi + bruit_phi;

        qubit_correction.set_theta(theta_corrige);
        qubit_correction.set_phi(phi_corrige);
        qubit_correction.normalize();

        complexe alpha_corr = qubit_correction.get_alpha();
        complexe beta_corr = qubit_correction.get_beta();

        rk4_step(alpha_corr, beta_corr, t, dt, omega, omega_0, omega_1);

        qubit_correction.set_alpha(alpha_corr);
        qubit_correction.set_beta(beta_corr);
        qubit_correction.normalize();

        ancien_bruit_theta = bruit_theta;
        ancien_bruit_phi = bruit_phi;

        t += dt;
    }
}

fichier.close();
}

void simulation_multiple_noise_levels_repete(qubit q_init, double omega, double omega_0, double omega_1,
    double dt, double T_init, double T_max, int n, int repetitions) {
    std::vector<double> noise_levels = {0.01, 0.02, 0.03, 0.04};
    std::ofstream fichier("simulation_results_repeated.csv");

    if (!fichier.is_open()) {
        std::cerr << "Erreur: Impossible d'ouvrir le fichier simulation_results.csv" << std::endl;
        return;
    }

    fichier << "Noise Level,Repetition,Time,abs_alpha_bruit,abs_beta_bruit,abs_alpha_corrige,abs_beta_corrige" << std::endl;

    for (double noise_level : noise_levels) {
        for (int rep = 0; rep < repetitions; ++rep) {
            qubit qubit_bruit = q_init;
            qubit qubit_correction = q_init;

            double t = T_init;
            double bruit_theta = 0.0;
            double bruit_phi = 0.0;
            double ancien_bruit_theta = 0.0, ancien_bruit_phi = 0.0;

            for (int i = 0; i < n - 1; i++) {
                fichier << noise_level << "," << rep << "," << t << ","
                        << norm(qubit_bruit.get_alpha()) << "," << norm(qubit_bruit.get_beta()) << ","
                        << norm(qubit_correction.get_alpha()) << "," << norm(qubit_correction.get_beta()) << std::endl;

                complexe alpha = qubit_bruit.get_alpha();
                complexe beta = qubit_bruit.get_beta();

                rk4_step(alpha, beta, t, dt, omega, omega_0, omega_1);

                qubit_bruit.set_alpha(alpha);
                qubit_bruit.set_beta(beta);
                qubit_bruit.normalize();

                appliquer_bruit_aleatoire(qubit_bruit, noise_level, &bruit_theta, &bruit_phi);

                double theta_corrige = qubit_correction.get_theta() - ancien_bruit_theta + bruit_theta;
                double phi_corrige = qubit_correction.get_phi() - ancien_bruit_phi + bruit_phi;

                qubit_correction.set_theta(theta_corrige);
                qubit_correction.set_phi(phi_corrige);
                qubit_correction.normalize();

                complexe alpha_corr = qubit_correction.get_alpha();
                complexe beta_corr = qubit_correction.get_beta();

                rk4_step(alpha_corr, beta_corr, t, dt, omega, omega_0, omega_1);

                qubit_correction.set_alpha(alpha_corr);
                qubit_correction.set_beta(beta_corr);
                qubit_correction.normalize();

                ancien_bruit_theta = bruit_theta;
                ancien_bruit_phi = bruit_phi;

                t += dt;
            }
        }
    }

    fichier.close();
}


/*
Preparation des qubits dans l'état |+>
*/





// fonction de préparation d'un qubit dans l'état souhaité
// très similaire à la simulation en champ oscillant
// retourne un qubit

qubit preparation_etat_plus(qubit q_init, double omega , double omega_0 , double omega_1, double dt, double T_init, double T_max, int n){
    // On laisse tourner la simulation pendant 3 ou 4 périodes, puis on coupe.
    // Cependant, au lieu de simplement stopper la simulation, on force omega_1 = 0 --> champs statique

    std::vector<double> t(n);
    std::vector<complexe> alpha(n), beta(n);
    // On a du gâchis d'allocation de mémoire, la majorité des cases mémoires de nos vecteurs ne vont pas être
    // utilisées

    // Conditions initiales
    t[0] = T_init ;
    alpha[0] = q_init.get_alpha();
    beta[0] = q_init.get_beta();

    // définition de n_min, nombre d'étape avant la sélection

    int n_min = std::round((T_max-T_init) / (2*dt)); // laisse tourner la moitié de la simulation

    int rang=0; // pour obtenir les dernière valeurs de alpha et beta de l'itération



    std::ofstream fichier("preparation_qubit.csv");
    fichier << "Time alpha beta abs_alpha2 abs_beta2" << std::endl;

        
        for (int i = 0; i < n - 1; i++) {
            fichier << t[i] << " " << alpha[i] << " " << beta[i] << " " 
                    << norm(alpha[i]) << " " << norm(beta[i]) << std::endl;

            if ((0.5 - 5e-3 <= norm(alpha[i]) && norm(alpha[i]) <= 0.5 + 5e-3) && i > n_min){ //la précision ici est arbitraire
                rang = i;
                //Au lieu de break, on stoppe les oscillations du champs
                omega_1 = 0;
            };
        
        

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



// Fonction de simulation avec champ magnétique oscillant et bruit, avec correction et arrêt du champ
qubit preparation_bruit_correction(qubit q_init, double omega, double omega_0, double omega_1, double dt, double T_init, double T_max, int n, double noise_level) {

    // Qubit bruité (expérience principale)
    qubit qubit_bruit = q_init;
    // Qubit corrigé
    qubit qubit_correction = q_init;

    double t = T_init;
    double bruit_theta = 0.0;
    double bruit_phi = 0.0;
    double ancien_bruit_theta = 0.0, ancien_bruit_phi = 0.0; // Bruit de l'itération précédente
    int rang;
    bool condition_prepa = false;

    std::ofstream fichier("oscillating_magnetic_field_with_noise_and_correction.csv");

    // Vérifier si le fichier est ouvert
    if (!fichier.is_open()) {
        std::cerr << "Erreur: Impossible d'ouvrir le fichier oscillating_magnetic_field_with_noise_and_correction.csv" << std::endl;
        return q_init;
    }

    // Écrire l'en-tête
    fichier << "Time abs_alpha_bruit abs_beta_bruit abs_alpha_corrige abs_beta_corrige bruit_alpha bruit_beta" << std::endl;

    // Définition de n_min, nombre d'étapes avant la sélection
    int n_min = std::round((T_max - T_init) / (2 * dt));

    // Boucle de simulation
    for (int i = 0; i < n - 1; i++) {
        // Écrire les données dans le fichier
        fichier << t << " " << norm(qubit_bruit.get_alpha()) << " " << norm(qubit_bruit.get_beta()) << " "
                << norm(qubit_correction.get_alpha()) << " " << norm(qubit_correction.get_beta()) << " "
                << bruit_theta << " " << bruit_phi << std::endl;

        // Intégration RK4
        complexe alpha = qubit_bruit.get_alpha();
        complexe beta = qubit_bruit.get_beta();
        rk4_step(alpha, beta, t, dt, omega, omega_0, omega_1);

        // Mise à jour du qubit bruité
        qubit_bruit.set_alpha(alpha);
        qubit_bruit.set_beta(beta);
        qubit_bruit.normalize();

        // Appliquer le bruit
        appliquer_bruit_aleatoire(qubit_bruit, noise_level, &bruit_theta, &bruit_phi);

        // Correction du bruit
        double theta_corrige = qubit_correction.get_theta() - ancien_bruit_theta + bruit_theta;
        double phi_corrige = qubit_correction.get_phi() - ancien_bruit_phi + bruit_phi;
        qubit_correction.set_theta(theta_corrige);
        qubit_correction.set_phi(phi_corrige);
        qubit_correction.normalize();

        complexe alpha_corr = qubit_correction.get_alpha();
        complexe beta_corr = qubit_correction.get_beta();
        rk4_step(alpha_corr, beta_corr, t, dt, omega, omega_0, omega_1);

        qubit_correction.set_alpha(alpha_corr);
        qubit_correction.set_beta(beta_corr);
        qubit_correction.normalize();

        // Vérifier si |alpha|^2 est proche de 0.5 après n_min étapes et stopper les oscillations du champ
        if ((0.5 - 5e-3 <= norm(alpha_corr) && norm(alpha_corr) <= 0.5 + 5e-3) && i > n_min && !condition_prepa) {
            rang = i;
            std::cout << "rang " << rang << std::endl << "temps :" << t << std::endl;
            omega_1 = 0;
            condition_prepa = true;
        }

        // Mettre à jour les anciens bruits pour la prochaine itération
        ancien_bruit_theta = bruit_theta;
        ancien_bruit_phi = bruit_phi;

        // Mettre à jour le temps
        t += dt;
    }

    // Fermer le fichier
    fichier.close();

    // Retourner le qubit préparé
    return qubit_correction;
}