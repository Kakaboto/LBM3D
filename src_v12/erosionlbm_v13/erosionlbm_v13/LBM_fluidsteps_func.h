#pragma once
#pragma once
#include "LBM_input.h"
#include "LBM_classes.h"

void IC(const momentum_direction& e, direction_density& f, direction_density& ftemp, Solid_list& solid_list);

void macrovariables(velocity& u, density& rho, Solid_list& solid_list, direction_density& f, momentum_direction& e);


void edf(Solid_list& solid_list, velocity& u, density& rho, EDF& feq, momentum_direction& e, int forcedirection);
double calcfeq(density& rho, momentum_direction& e, double ueq[3], int ix, int iy, int iz, int a, int forcedirection);
double * calcforce(double ueq[3], density& rho, int ix, int iy, int iz, int forcedirection);
double edfvecdot(momentum_direction& e, double ueq[3], int a);

void stream(Solid_list& solid_list, direction_density& f, direction_density& ftemp, momentum_direction& e);
void bouncebackBC(direction_density& ftemp, direction_density& f, momentum_direction& e, int ix, int iy, int iz, int a, int ashift);

void collision(Solid_list& solid_list, direction_density& f, direction_density& ftemp, EDF& feq, momentum_direction& e);