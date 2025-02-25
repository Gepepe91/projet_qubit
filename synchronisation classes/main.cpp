#include<iostream>
#include<cmath>
#include<fstream>
#include"matrice.h"
#include<complex>
#include"qubit.h"

using complexe = std::complex<double>;
using namespace std;

int main(){

    qubit q1(complexe(1,M_PI) , complexe(1,10)) , q2(complexe(1,-1) , complexe(-1,1));

    matrice m;

    m = q1 ^ q2;

    m.display();

    return 0;
}