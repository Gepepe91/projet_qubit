#include<iostream>
#include<cmath>
#include"complexe.h"

//private


void complexe::i_fois(){
    // x --> -y
    // y --> x
    double x = part_reelle ; double y = part_imaginaire;
    part_reelle = - y;
    part_imaginaire = x;
}

void complexe::synchr_pol_to_cart(){
    //mise à jour de la forme polaire à partir de la forme cartésienne
    part_reelle = attr_norme * cos(attr_argument);
    part_imaginaire = attr_norme * sin(attr_argument);
}

void complexe::synchr_cart_to_pol(){
    //mise à jour de la forme cartésienne à partir de la forme polaire
    attr_norme = sqrt(part_reelle*part_reelle + part_imaginaire*part_imaginaire);
    attr_argument = atan2(part_imaginaire , part_reelle);
}

//public

complexe::complexe(){
    part_reelle=0;
    part_imaginaire=0;
}

complexe::complexe(double x_ , double y_) : part_reelle(x_) , part_imaginaire(y_) {}



double complexe::get_part_reelle(){
    synchr_pol_to_cart();
    return part_reelle;
}

double complexe::get_part_imaginaire(){
    synchr_pol_to_cart();
    return part_imaginaire;
}

void complexe::initialiser_avec_norme_argument(double r, double theta){
    attr_norme = r;
    attr_argument = theta;
    synchr_pol_to_cart();
}

void complexe::initialiser(double x , double y){
    part_reelle = x;
    part_imaginaire = y;
    synchr_cart_to_pol();
}

void complexe::afficher(){ 
    std::cout << part_reelle << " + " << "i(" << part_imaginaire << ')' ;//<< std::endl;
}

double complexe::get_norm(){
    synchr_cart_to_pol();
    return attr_norme;
}

double complexe::get_phase(){
    synchr_cart_to_pol();
    return attr_argument;
}
