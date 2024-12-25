#include <iostream>
#include <complex>
#include <vector>
#include <fstream>

using namespace std;

// Définir un type pour un vecteur à 2 dimensions
using Vector2 = complex<double>[2]; //dans la classe qubit, cette forme plus compact peut être bien

// Fonction pour multiplier une matrice 2x2 par un vecteur 2D
void matVecMul(const complex<double> H[2][2], const complex<double> psi[2], complex<double> result[2]) {
    result[0] = H[0][0] * psi[0] + H[0][1] * psi[1];
    result[1] = H[1][0] * psi[0] + H[1][1] * psi[1];
}

// Fonction pour calculer f(t, psi) = -i * H * psi
void f(const complex<double> H[2][2], const complex<double> psi[2], complex<double> result[2]) {
    matVecMul(H, psi, result);
    result[0] *= complex<double>(0, -1); // Multiplier par -i
    result[1] *= complex<double>(0, -1); // Multiplier par -i
}

// Méthode de Runge-Kutta 4
void rungeKutta4(const complex<double> H[2][2], complex<double> psi[2], double dt) {
    complex<double> k1[2], k2[2], k3[2], k4[2], temp[2];

    // k1
    f(H, psi, k1);

    // k2
    temp[0] = psi[0] + 0.5 * dt * k1[0];
    temp[1] = psi[1] + 0.5 * dt * k1[1];
    f(H, temp, k2);

    // k3
    temp[0] = psi[0] + 0.5 * dt * k2[0];
    temp[1] = psi[1] + 0.5 * dt * k2[1];
    f(H, temp, k3);

    // k4
    temp[0] = psi[0] + dt * k3[0];
    temp[1] = psi[1] + dt * k3[1];
    f(H, temp, k4);

    // Mise à jour de psi
    psi[0] += (dt / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]);
    psi[1] += (dt / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]);
}

int main() {
    // Hamiltonien H
    complex<double> H[2][2] = {
        {complex<double>(0, 0), complex<double>(-1, 0)},
        {complex<double>(-1, 0), complex<double>(0, 0)}
    };

    // État initial |psi(0)>
    complex<double> psi[2] = {1.0, 0.0}; // Par exemple |0>

    // Paramètres de la simulation
    double t_max = 10.0;
    double dt = 0.01;

    //création et ouverture du fichier
    ofstream fich;
    fich.open("test.csv");

    // Simulation
    cout << "Time\tReal(psi0)\tImag(psi0)\tReal(psi1)\tImag(psi1)" << endl;
    for (double t = 0; t <= t_max; t += dt) {
        //affichage dans le terminal (optionnel)
        cout << t << "\t" 
             << psi[0].real() << "\t" << psi[0].imag() << "\t" 
             << psi[1].real() << "\t" << psi[1].imag() << endl;
        // écriture dans le fichier
        fich << t << "\t" 
             << psi[0].real() << "\t" << psi[0].imag() << "\t" 
             << psi[1].real() << "\t" << psi[1].imag() << endl;
        rungeKutta4(H, psi, dt); // Mise à jour avec RK4
    }

    //fermeture du fichier
    fich.close();

    return 0;
}
