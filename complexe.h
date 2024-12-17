#ifndef COMPLEXES_H
#define COMPLEXES_H

class complexe
{

public:

complexe();
complexe(double x_ , double y_);

double get_part_reelle();
double get_part_imaginaire();

void initialiser(double x , double y);
void initialiser_avec_norme_argument(double r , double theta);
void afficher();
double get_norm();
double get_phase();






private:

double part_reelle , part_imaginaire ;
double attr_norme , attr_argument; // attributs privés, mais on peut les obtenir avec une méthode

void synchr_cart_to_pol();
void synchr_pol_to_cart();

void i_fois();

};

#endif