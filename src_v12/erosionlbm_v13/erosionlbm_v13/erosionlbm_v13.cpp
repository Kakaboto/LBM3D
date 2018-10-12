// erosionlbm.cpp : Defines the entry point for the console application.
//
#pragma once

#include "LBM_classes.h"
#include "LBM_input.h"
#include "LBM_force_erosion_func.h"
#include "LBM_fluidsteps_func.h"
#include "LBM_utility_func.h"

//... should prob försöka lägga alla vektorer i en å samma klass med massa olika () operatorer bara.




int main()
{
	//===========================================================================================================
	// Initializing variables
	FILE * velfile = fopen("velocity_v9.txt", "w");
	FILE * densfile = fopen("density_v9.txt", "w");
	FILE * gridfile = fopen("grid_v9.txt", "w");
	FILE * parfile = fopen("parameters_v9.txt", "w");
	FILE * solfile = fopen("solid_list_v9.txt", "w");
	FILE * reyfile = fopen("rey_umax_v9.txt", "w");
	FILE * stressfile = fopen("shearstress_v9.txt", "w");
	FILE * sttensfile = fopen("stresstensor_v9.txt", "w");
	FILE * forcefile = fopen("force_v9.txt", "w");
	FILE * nhatfile = fopen("nhat_v9.txt", "w");
	FILE * torfile = fopen("torfile_v9.txt", "w");
	FILE * erodefile = fopen("erodefile_v9.txt", "w");
	FILE * errorfile = fopen("errorfile_v9.txt", "w");
	FILE * eronumbfile = fopen("erosionnumber_v9.txt", "w"); // prints out erosionnumber squared.
	FILE * dmfile = fopen("dm_v9.txt", "w");
	int printi;
	double umax = 0;
	double uav = 0;

	Grid grid;
	momentum_direction e;
	EDF feq;
	velocity u;
	density rho; // IC: rho = 1 on all fluid nodes.
	direction_density f;
	direction_density ftemp;
	Stresstensor stresstensor; //nearest neighbour to surface solid nodes.
	Normalvector nhat;
	Wall_force tau_stress;
	vector3Ncubed F_D;
	vector3Ncubed F_sum;
	vector3Ncubed torque;
	density masschange;
	vectorNcubed ero_reso_check; //erosion resolution check. It's used to check if we erode faster than Delta_T or not.
	vectorNcubed F_vdw;
	masschange.clear();
	vector<int> BCtype = readBC(x0BC, xNBC, y0BC, yNBC, z0BC, zNBC);
	//==========================================================================================================
	// Choosing object to put in flow and applying IC.
	grid.printgrid(gridfile);
	cout << "xmax = " << (Nx - 1)*latspace << ", ymax = " << (Ny - 1)*latspace << ", zmax = " << (Nz - 1)*latspace << "\n";

	cout << "Choose object. 1: sphere. 2: cylinder pipe. 3: square pipe. 4: cylinder with sphere inside. 5: cylinder. 6: open box. 7: triangle-like. \n";
	//cin >> obchoice;
	//for (double rotation = -0.25; rotation <  0.26; rotation += 0.05) {
	Solid_list solid_list(obchoice, grid, pi * 0, e, parfile);
	updatePBC_solid(solid_list);
	solid_list.printsolid_list(solfile);

	// IC for rho (rho = 1 everywhere).
	//	for (int ix = 0; ix < Nx; ix++) {
	//		for (int iy = 0; iy < Ny; iy++) {
	//			for (int iz = 0; iz < Nz; iz++) {
	//				rho(ix, iy, iz) = 1.;
	//			}
	//		}
	//	}
	// Kolla att figuren som printades ut i solid_list roterar korrekt. Sen kör en lång simulering med olika vinklar.
	cout << "Grid and object created. Running simulation for t = \n";
	//cin >> tend;
	IC(e, f, ftemp, solid_list);
	updateBC(f, -1, Bvel, rho, e, u, BCtype);
	updateBC(ftemp, -1, Bvel, rho, e, u, BCtype);
	printi = 0;
	//==========================================================================================================
	// Main program.
	int i_er = 1;
	int i_Fvdw = 1;

	for (int t = 0; t < tend; t++) {
		cout << t << " ";
		stream(solid_list, f, ftemp, e);
		updateBC(ftemp, t, Bvel, rho, e, u, BCtype);
		macrovariables(u, rho, solid_list, ftemp, e);
		//		umax = find_umax(u);
		//		uav = find_uav(u);
		//		fprintf(reyfile, "%e %e %e\n", umax, Bvel[0] * 16. / mu, uav);
		if (t == Delta_T * printi) {
			solid_list.printsolid_list(solfile);
			printstuff(velfile, densfile, parfile, reyfile, stressfile, forcefile, nhatfile, sttensfile, torfile, erodefile, eronumbfile, dmfile, t, u, rho, tau_stress, F_D, nhat, stresstensor, torque, masschange, F_vdw, solid_list, masschange);
			printi++;
		}
		edf(solid_list, u, rho, feq, e, edfforcedir);
		collision(solid_list, f, ftemp, feq, e);
		updateBC(f, t, Bvel, rho, e, u, BCtype);
		computestress(e, ftemp, f, feq, solid_list, stresstensor, nhat, tau_stress, F_sum, masschange, F_vdw, i_er, i_Fvdw, rho);
//		if (i_Fvdw == i_er)
//			i_Fvdw++;
		computetorque(solid_list, tau_stress, torque);
//		if (t == i_er*Delta_T) {
		erosion(solid_list, e, F_sum, rho, f, edfforcedir, solfile, masschange, nhat, ero_reso_check, F_vdw, errorfile);
//			i_er++;
//		}
	}
	cout << "\n Done! \n";
	solid_list.clear_list();
	//	};
	//==========================================================================================================


	//	getchar();
	fclose(velfile);
	fclose(densfile);
	fclose(gridfile);
	fclose(parfile);
	fclose(solfile);
	fclose(forcefile);
	fclose(stressfile);
	fclose(torfile);

	return 0;
}









/*
void updateBC_v2(direction_density& f, int t, double Bvel[3], density& rho, momentum_direction& e, velocity& u) {
int ix = 0;
int iy = 0;
int iz = 0;

for (iz = -1; iz < Nz + 1; iz++) {
for (iy = -1; iy < Ny + 1; iy++) {
for (ix = -1; ix < Nx + 1; ix++) {
if (isnodeBC(ix, iy, iz) == 1) {

}
}
}
}
}
int isnodeBC(int ix, int iy, int iz) {
if (ix == -1 || ix == Nx || iy == -1 || iy == Ny || iz == -1 || iz == Nz)
return 1;
else
return 0;
}
*/
