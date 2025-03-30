#include "fonctions.h"
#include "constants.h"
#include "qubit.h"
#include <complex>
#include <cmath>

using complexe = std::complex<double>;

qubit compute_derivative(qubit q, matrice H) {
    // Appliquer la matrice H au qubit
    qubit new_q = H*q;
    // multiplication par -i / hbar
    new_q = new_q * (complexe (0,-1) * (1/hbar));
    return new_q;
}

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

void simulation_static(qubit q , matrice H , double dt , double T_init, double T_max, int n){
    //Ouverture du fichier
    std::ofstream fichier;
    fichier.open("data.csv");
    double t = T_init;

    //boucle for, simulation à proprement parlé
    fichier << "Time theta phi abs_alpha2 abs_beta2" << std::endl;
    for (int i = 0; i < n - 1; i++) {
        fichier << t  << " " << q.get_theta() << " " << q.get_phi() << " " << q.get_abs_alpha2() << " " << q.get_abs_beta2() << std::endl;
        rk4(q, H, dt); // Mise à jour de |psi> = notre qubit avec RK4
        q.normalize(); // renormalisation, sinon le qubit se dénormalise au cours du temps (erreur numérique)
        t = t + dt;
    };
    //Fermeture du fichier
    fichier.close();
}









// Définition des équations différentielles couplées
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
    fichier << "Time abs_alpha2 abs_beta2" << std::endl;

    for (int i = 0; i < n - 1; i++) {
        fichier << t[i] << " " << norm(alpha[i]) << " " << norm(beta[i]) << std::endl;

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



static std::random_device rd;
static std::mt19937 gen(rd());


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
    std::cout << "Fichier data_with_noise.csv créé avec succès." << std::endl; // Message de débogage
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