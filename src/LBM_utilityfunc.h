#pragma once
#include "LBM_classes.h"
#include "../LBM_globvar.h"

dvec crossproduct(Wall_force& tau_stress, double r[3], int ix, int iy, int iz);

vector<int> readBC(string x0BC, string xNBC, string y0BC, string yNBC, string z0BC, string zNBC);

void printstuff(FILE * velfile, FILE * densfile, FILE * parfile, FILE * reyfile, FILE * stressfile, FILE * forcefile, FILE * nhatfile, FILE * sttensfile, FILE * torfile, FILE * erodefile, int Renum, velocity& u, density& rho, Wall_force& tau_stress, vector3Ncubed& F_D, Normalvector& nhat, Stresstensor& stresstensor, vector3Ncubed& torque, density& erodelist);

void updatePBC(direction_density& f);
void updateBC(direction_density& f, int t, double Bvel[3], density& rho, momentum_direction& e, velocity& u,  vector<int>& BCtype);
void updatePBC_solid(Solid_list& solid_list);
