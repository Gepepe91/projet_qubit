#include<fstream>

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









qubit preparation_etat_plus(qubit q_init, double omega , double omega_0 , double omega_1, double dt, double T_init, double T_max, int n) {
    std::vector<double> t(n);
    std::vector<complexe> alpha(n), beta(n);
    // On a du gâchis d'allocation de mémoire, la majorité des cases mémoires de nos vecteurs ne vont pas être
    // utilisées

    // Conditions initiales
    t[0] = T_init ;
    alpha[0] = q_init.get_alpha();
    beta[0] = q_init.get_beta();



    std::ofstream fichier("preparation_qubit.csv");
    fichier << "Time alpha beta abs_alpha2 abs_beta2" << std::endl;

        
        for (int i = 0; i < n - 1; i++) {
            fichier << t[i] << " " << alpha[i] << " " << beta[i] << " " 
                    << norm(alpha[i]) << " " << norm(beta[i]) << std::endl;

            if (0.5 - 1e-5 <= norm(alpha[i]) && norm(alpha[i]) <= 0.5 + 1e-5) break;

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
}
