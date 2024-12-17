#include"complexe.h" //techniquement inutile, simplement pour calmer VSCODE qui panique

//Entre complexe

complexe operator+(complexe z1 , complexe z2){
    complexe z3;
    z3.initialiser(z1.get_part_reelle() + z2.get_part_reelle() , z1.get_part_imaginaire() + z2.get_part_imaginaire());
    return z3;
}

complexe operator-(complexe z){ // version unaire, inversion de signe
    complexe _z;
    _z.initialiser(-z.get_part_reelle() , -z.get_part_imaginaire());
    return _z;
}

complexe operator-(complexe z1 , complexe z2){ // version binaire, operatrion entre deux complexe
    complexe z3;
    z3.initialiser(z1.get_part_reelle() - z2.get_part_reelle() , z1.get_part_imaginaire() - z2.get_part_imaginaire());
    return z3;
}

complexe operator*(complexe z1 , complexe z2){
    complexe z3;
    z3.initialiser_avec_norme_argument(z1.get_norm() * z2.get_norm() , z1.get_phase() + z2.get_phase());
    return z3;
}

complexe operator/(complexe z1 , complexe z2){
    complexe z3;
    z3.initialiser_avec_norme_argument(z1.get_norm() / z2.get_norm() , z1.get_phase() - z2.get_phase());
    return z3;
}

//complexe et réels
//Les opérateurs d'addition/soustraction prennent le complexe à gauche et le réel à droite
//Les opérateurs de multiplication/division prennent le réel à gauche et le complexe à droite

complexe operator+(complexe z , double x){
    complexe zz;
    zz.initialiser(z.get_part_reelle() + x , z.get_part_imaginaire());
    return zz;
}

complexe operator-(complexe z , double x){
    complexe zz;
    zz.initialiser(z.get_part_reelle() - x , z.get_part_imaginaire());
    return zz;
}

complexe operator*(double x , complexe z){
    complexe zz;
    zz.initialiser_avec_norme_argument(x*z.get_norm() , z.get_phase());
    return zz;
}

complexe operator/(double x , complexe z){
    complexe zz;
    zz.initialiser_avec_norme_argument(x/z.get_norm() , -z.get_phase());
    return zz;
}