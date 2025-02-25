#include<iostream>
#include<cmath>
#include<fstream>
#include"matrice.h"
#include<complex>
#include"qubit.h"

using complexe = std::complex<double>;
using namespace std;

int main(){

    matrice m;

    m.dephasage(M_PI/2);

    m.display();

    return 0;
}