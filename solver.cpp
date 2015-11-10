#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

class solver{
	/************** Paramètres ****************/
private:
	//dimension
	int n;
	//taille de la grille
	double L;
	//nombre de cellules 
	int N;
	//dimension d'un voxel
	double D = L / N;
	//intervalle de temps
	double dt;
	//grilles sources : quel format ?
	Mat S0, S1;
	//grilles sur lesquelles on travaille
	Mat U1, U0;
	//viscosité
	double visc;
	//coefficient de diffusion
	double k;
	//coefficient de dissipation 
	double a;

};