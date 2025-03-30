#ifndef TENS_2Q_H
#define TENS_2Q_H

/*
Dans l'idée, on veut créer une classe fille de qubit : la classe tens_2q.
C'est classe englobe les combinaisons de 2 qubits, ainsi que 2 qubits intriqués.
Elle doit hériter des attributs alpha et beta, mais elle ajoute gamma et delta.

Par contre, la surcharge d'opérateur avec des matrices 2x2 n'a plus de sens.
Il faut en redéfinir (en s'inspirant grandement) pour des matrices H1 tens H2 de dimensions 4x4
*/

#include "qubit.h"

class tens_2q : public qubit { //|00>, |01>, |10>, |11> ; choix arbitraire
protected:

    std::complex<double> gamma, delta; // Amplitudes pour |10> et |11>

public:
    // Constructeurs
    tens_2q();
    tens_2q(std::complex<double> alpha_, std::complex<double> beta_,
                    std::complex<double> gamma_, std::complex<double> delta_);

    // Getters
    std::complex<double> get_gamma() const;
    std::complex<double> get_delta() const;
    double get_abs_gamma2();
    double get_abs_delta2();

    // Setters
    void set_gamma(std::complex<double> gamma_);
    void set_delta(std::complex<double> delta_);

    // Méthode de normalisation spécifique aux deux qubits intriqués
    void normalize_entrelace();

    // Affichage des amplitudes des deux qubits
    void display() const;
};

#endif
