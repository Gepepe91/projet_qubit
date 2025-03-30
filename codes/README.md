# Simulation de Qubit en C++

Un programme écrit en C++ permettant de simuler le comportement de systèmes de qubits dans des champs oscillants.

## Description

Dans ce programme, il est possible de simuler le comportement d'un qubit plongé dans un champ statique et dans un champ dynamique.

## Installation

Se placer dans le répertoire, puis entrer :
"make"


## Utilisation

./programme



## Structure du projet
```
mon-projet/
│── src/               # Fichiers sources C++
│   ├── qubit.cpp
│   ├── fonctions.cpp
│   └── ...
│── include/           # Headers des classes
│   ├── qubit.h
│   ├── fonctions.h
│   └── ...
│── main.cpp           # Point d'entrée du programme
│── Makefile           # Script de compilation
│── README.md          # Documentation du projet
│── visualisation.ipynb# Notebook Python pour obtenir les figures
│── data/              # Dossier contenant les données en cas de soucis avec le code
    ├── simulation_statique.csv
    ├── simulation_dynamique.csv
    ├── ...
```

## Fichiers de sorties

Nous avons d'autres fichiers de données créés par notre programme. Ici, nous résumons les principaux :

simulation_statique.csv     #données de la simulation d'un qubit dans un champ statique
simulation_dynamique.csv    #données de la simulation d'un qubit ndans un champ dynamique
preparation_qubit.csv       # données de la préparation d'un qubit dans l'état |+>
preparation_qubit_champ_dynamique_bruit_correction


##  Auteurs

Gabriel Pus Perchaud
Hassan Id Taleb
Florien Giletta


## Licence
Ce projet n'est pas protégé 🥲
