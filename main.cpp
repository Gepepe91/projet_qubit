#include<iostream>
#include<cmath>
#include<sstream>
#include"complexe.h"
#include"surchage_complexe.h"
#include"qubit.h"

using namespace std;

int main(){

    complexe alpha(1,1) , beta(-10,1);

    qubit q(alpha , beta);


    q.display();
    
return 0;
}