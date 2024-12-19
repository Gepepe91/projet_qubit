#include<iostream>
#include<cmath>
//#include<Python.h>
#include<sstream>
#include"complexe.h"
#include"surchage_complexe.h"
#include"qubit.h"

using namespace std;

int main(){

    complexe alpha(1,1) , beta(-10,1);

    qubit q(alpha , beta);


    q.display();

    ostringstream ostpy;

    ostpy
    << "print('fuck')\n"
    << "import numpy as np\n"
    << "print(np.arange(10))\n"
    << "print(alpha.get_norm())"
        ;

    Py_Initialize();
    PyRun_SimpleString(ostpy.str().c_str());
    Py_Finalize();
    
return 0;
}