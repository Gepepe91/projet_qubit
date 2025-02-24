#include "Matrice.h"

using namespace std;

int main(){
    Matrice hadamard(2) ;
    Matrice PauliZ(2) ;
    hadamard.gateHADAMARD() ;
    hadamard.affichermatrice() ;
    PauliZ.PauliZ() ;
    PauliZ.affichermatrice() ;

    return 0 ;

}