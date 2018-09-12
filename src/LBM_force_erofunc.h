#pragma once
#include "LBM_classes.h"
#include "../LBM_globvar.h"
#include "LBM_utilityfunc.h"

void computestress(momentum_direction& e, direction_density& ftemp, direction_density& f, EDF& feq, Solid_list& solid_list, Stresstensor& stresstensor, Normalvector& nhat, Wall_force& tau_stress, vector3Ncubed& F_D);

void computetorque(Solid_list& solid_list, Wall_force& tau_stress, vector3Ncubed& torque);

void erosion(Solid_list& solid_list, momentum_direction& e, vector3Ncubed& F_sum, density& rho, direction_density& f, int forcedirection, FILE * solfile, density& erodelist, Normalvector& nhat);
