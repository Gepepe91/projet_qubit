#include "tens_Nq.h"
#include <iostream>
#include <cmath>

// Constructeur
tens_Nq::tens_Nq(int n_) : n(n_) {
    int dim = 1 << n; // 2^n
    amplitudes.resize(dim, std::complex<double>(0, 0));
    amplitudes[0] = 1.0; // État initial par défaut |00...0>
}

// Normalisation
void tens_Nq::normalize() {
    double norm = 0.0;
    for (const auto& amp : amplitudes) {
        norm += std::norm(amp);
    }
    norm = std::sqrt(norm);
    for (auto& amp : amplitudes) {
        amp /= norm;
    }
}

// Affichage de l'état
void tens_Nq::display() const {
    std::cout << "État du système à " << n << " qubits :" << std::endl;
    for (size_t i = 0; i < amplitudes.size(); ++i) {
        std::cout << "|";
        for (int j = n - 1; j >= 0; --j) {
            std::cout << ((i >> j) & 1);
        }
        std::cout << "> : " << amplitudes[i] << std::endl;
    }
}

// Produit scalaire
std::complex<double> tens_Nq::operator|(const qubit_system& qs) const {
    const tens_Nq* q = dynamic_cast<const tens_Nq*>(&qs);
    if (q && q->n == n) {
        std::complex<double> result = 0;
        for (size_t i = 0; i < amplitudes.size(); ++i) {
            result += std::conj(amplitudes[i]) * q->amplitudes[i];
        }
        return result;
    } else {
        throw std::runtime_error("Produit scalaire non défini pour ce type.");
    }
}


// On admet, c'est du chatGPT parce que c'est prise de tête
matrice tens_Nq::operator&(const qubit_system& qs) const {
    // Convertir le "qs" en un tens_Nq pour effectuer l'opération
    const tens_Nq* q = dynamic_cast<const tens_Nq*>(&qs);
    if (!q) {
        throw std::runtime_error("Produit tensoriel non défini pour ce type.");
    }

    // Préparez un état tensoriel qui représente le produit des deux systèmes
    std::vector<std::complex<double>> new_amplitudes(1 << (this->get_n() + q->get_n()));  // Taille = 2^(N+M)

    // Parcours de chaque amplitude du premier état et du deuxième état
    for (size_t i = 0; i < (1 << this->get_n()); ++i) {
        for (size_t j = 0; j < (1 << q->get_n()); ++j) {
            new_amplitudes[(i << q->get_n()) + j] = this->amplitudes[i] * q->amplitudes[j];
        }
    }

    // Créer et retourner le résultat sous forme de matrice
    // La matrice doit avoir 2^(n+m) lignes et 1 colonne pour correspondre à un vecteur
    matrice result(1 << (this->get_n() + q->get_n()), 1);
    for (size_t i = 0; i < new_amplitudes.size(); ++i) {
        result.setElement(i, 0, new_amplitudes[i].real());
    }

    return result;
}
