#include<iostream>
#include<cmath>
#include<sstream>
#include<complex>
#include"qubit.h"
#include"matrice.h"

using namespace std;

int main(){

    complex<double> a(4.,4) , b(0,0) , c(-5,5) , d(1,1);

    matrice m(a , b , c , d);

    m.display();

    m.get_element(0,1).afficher();

    cout << endl <<  endl;

    m.set_as_Pauli_z();
    m.display();

    return 0;
}