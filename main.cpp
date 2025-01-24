#include<iostream>
#include<cmath>
#include<sstream>
#include<complex>
#include"qubit.h"
#include"matrice.h"

using complexe = std::complex<double>;
using namespace std;

int main(){

    complexe a(4.,4) , b(0,0) , c(-5,5) , d(1,1);

    matrice m(a , b , c , d);

    m.display();

    cout << endl;

    cout << m.get_element(0,1);

    cout << endl <<  endl;

    m.set_as_Pauli_z();
    m.display();

    return 0;
}
