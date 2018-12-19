
#include "LBM_utility_func.h"

void updateBC(direction_density& f, int t, double Bvel[3], density& rho, momentum_direction& e, velocity& u, vector<int> BCtype) {
	int a = 0;
	int ix = 0;
	int iy = 0;
	int iz = 0;
	// 1 - periodic BC. 2 - velocity BC.
	for (a = 0; a < 27; a++) {
		if (BCtype[5] == 1) {	//periodic BC
								//----------------------------------------------------------------------
								// upper boundary layer 
								// main surface
			iz = 0;
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					f.element[f.index(ix, iy, iz + Nz, a)] = f(ix, iy, iz, a);
				}
			}
			// edges
			// copying to left upper edge
			//			for (iy = 0; iy < Ny; iy++)
			//				f.element[f.index(-1, iy, iz + Nz, a)] = f(Nx - 1, iy, iz, a);
			// copying to right upper edge
			//			for (iy = 0; iy < Ny; iy++)
			//				f.element[f.index(Nx, iy, iz + Nz, a)] = f(0, iy, iz, a);

			// copying to forward upper edge
			for (ix = 0; ix < Nx; ix++)
				f.element[f.index(ix, -1, iz + Nz, a)] = f(ix, Ny - 1, iz, a);

			// copying to back upper edge
			for (ix = 0; ix < Nx; ix++)
				f.element[f.index(ix, Ny, iz + Nz, a)] = f(ix, 0, iz, a);

			// corners
			//			f.element[f.index(-1, -1, Nz, a)] = f(Nx - 1, Ny - 1, 0, a);		// upper left-forward
			//			f.element[f.index(Nx, -1, Nz, a)] = f(0, Ny - 1, 0, a);			// upper right-forward
			//			f.element[f.index(Nx, Ny, Nz, a)] = f(0, 0, 0, a);	// upper right-back
			//			f.element[f.index(-1, Ny, Nz, a)] = f(Nx - 1, 0, 0, a);		// upper left-back
		}
		if (BCtype[5] == 2 && t == -1) {	// velocity BC, outlet
											//----------------------------------------------------------------------
											// upper boundary layer
											// main surface
			iz = 0;
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					f.element[f.index(ix, iy, iz + Nz, a)] = f.element[f.index(ix, iy, iz + Nz - 1, a)];
				}
			}
			// edges
			// left upper edge
			for (iy = 0; iy < Ny; iy++)
				f.element[f.index(-1, iy, iz + Nz, a)] = f.element[f.index(0, iy, iz + Nz - 1, a)];

			// right upper edge
			for (iy = 0; iy < Ny; iy++)
				f.element[f.index(Nx, iy, iz + Nz, a)] = f.element[f.index(Nx - 1, iy, iz + Nz - 1, a)];

			// forward upper edge
			for (ix = 0; ix < Nx; ix++)
				f.element[f.index(ix, -1, iz + Nz, a)] = f.element[f.index(ix, 0, iz + Nz - 1, a)];

			// back upper edge
			for (ix = 0; ix < Nx; ix++)
				f.element[f.index(ix, Ny, iz + Nz, a)] = f.element[f.index(ix, Ny - 1, iz + Nz - 1, a)];

			// corners
			f.element[f.index(-1, -1, Nz, a)] = f.element[f.index(0, 0, Nz - 1, a)];			// upper left-forward
			f.element[f.index(Nx, -1, Nz, a)] = f.element[f.index(Nx - 1, 0, Nz - 1, a)];		// upper right-forward
			f.element[f.index(Nx, Ny, Nz, a)] = f.element[f.index(Nx - 1, Ny - 1, Nz - 1, a)];	// upper right-back
			f.element[f.index(-1, Ny, Nz, a)] = f.element[f.index(0, Ny - 1, Nz - 1, a)];		// upper left-back
		}
		//----------------------------------------------------------------------
		// lower boundary layer (copiyng from real upper layer)
		// main surface
		if (BCtype[4] == 1) {	// Periodic BC
			iz = Nz - 1;
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					f.element[f.index(ix, iy, -1, a)] = f(ix, iy, iz, a);
				}
			}
			// edges
			// copying to left lower edge
			//			for (iy = 0; iy < Ny; iy++)
			//				f.element[f.index(-1, iy, -1, a)] = f(Nx - 1, iy, iz, a);

			// copying to right lower edge
			//			for (iy = 0; iy < Ny; iy++)
			//				f.element[f.index(Nx, iy, -1, a)] = f(0, iy, iz, a);

			// copying to forward lower edge
			for (ix = 0; ix < Nx; ix++)
				f.element[f.index(ix, -1, -1, a)] = f(ix, Ny - 1, iz, a);

			// copying to back lower edge
			for (ix = 0; ix < Nx; ix++)
				f.element[f.index(ix, Ny, -1, a)] = f(ix, 0, iz, a);

			// corners
			//			f.element[f.index(-1, -1, -1, a)] = f(Nx - 1, Ny - 1, Nz - 1, a);		// bottom left-forward
			//			f.element[f.index(Nx, -1, -1, a)] = f(0, Ny - 1, Nz - 1, a);		// bottom right-forward
			//			f.element[f.index(Nx, Ny, -1, a)] = f(0, 0, Nz - 1, a);			// bottom right-back
			//			f.element[f.index(-1, Ny, -1, a)] = f(Nx - 1, 0, Nz - 1, a);		// bottom left-back
		}
		if (BCtype[4] == 2 && t == -1) {	// Velocity BC, inlet
			iz = Nz - 1;
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					f.element[f.index(ix, iy, -1, a)] = calcfeq(rho, e, Bvel, 1, 1, 1, a, 0);
				}
			}
			// edges

			// copying to right lower edge
			for (iy = 0; iy < Ny; iy++)
				f.element[f.index(Nx, iy, -1, a)] = calcfeq(rho, e, Bvel, 1, 1, 1, a, 0);

			// copying to forward lower edge
			for (ix = 0; ix < Nx; ix++)
				f.element[f.index(ix, -1, -1, a)] = calcfeq(rho, e, Bvel, 1, 1, 1, a, 0);

			// copying to back lower edge
			for (ix = 0; ix < Nx; ix++)
				f.element[f.index(ix, Ny, -1, a)] = calcfeq(rho, e, Bvel, 1, 1, 1, a, 0);

			// corners
			f.element[f.index(-1, -1, -1, a)] = calcfeq(rho, e, Bvel, 1, 1, 1, a, 0);		// bottom left-forward
			f.element[f.index(Nx, -1, -1, a)] = calcfeq(rho, e, Bvel, 1, 1, 1, a, 0);		// bottom right-forward
			f.element[f.index(Nx, Ny, -1, a)] = calcfeq(rho, e, Bvel, 1, 1, 1, a, 0);			// bottom right-back
			f.element[f.index(-1, Ny, -1, a)] = calcfeq(rho, e, Bvel, 1, 1, 1, a, 0);		// bottom left-back
		}
		//-----------------------------------------------------------------------------------------------------------
		// set right side ghost points
		if (BCtype[1] == 1) {	// Periodic BC
			ix = Nx;
			for (iz = 0; iz < Nz; iz++) {
				f(Nx, -1, iz, a) = f(0, Ny - 1, iz, a);			// copy to right-forward corner side
				f(Nx, Ny, iz, a) = f(0, 0, iz, a);				// copy to right-back corner side
				for (iy = 0; iy < Ny; iy++) {
					f.element[f.index(ix, iy, iz, a)] = f(0, iy, iz, a);
				}
			}
		}
		if (BCtype[1] == 2 && t == -1) {	// Velocity BC, outlet
			ix = Nx;
			for (iz = 0; iz < Nz; iz++) {
				f(Nx, -1, iz, a) = f(Nx - 1, 0, iz, a);			// copy to right-forward edge
				f(Nx, Ny, iz, a) = f(Nx - 1, Ny - 1, iz, a);	// copy to right-back edge
				for (iy = 0; iy < Ny; iy++) {
					f.element[f.index(ix, iy, iz, a)] = f.element[f.index(ix - 1, iy, iz, a)];	// Copy f:s in outlet to ghost points from real points.
				}
			}

			// copying to right upper edge
			for (iy = 0; iy < Ny; iy++)
				f.element[f.index(Nx, iy, Nz, a)] = f(Nx - 1, iy, Nz - 1, a);
			// copying to right lower edge
			for (iy = 0; iy < Ny; iy++)
				f.element[f.index(Nx, iy, -1, a)] = f(Nx - 1, iy, 0, a);

			f.element[f.index(Nx, -1, Nz, a)] = f(Nx - 1, 0, Nz - 1, a);			// upper right-forward
			f.element[f.index(Nx, Ny, Nz, a)] = f(Nx - 1, Ny - 1, Nz - 1, a);		// upper right-back
			f.element[f.index(Nx, -1, -1, a)] = f(Nx - 1, 0, 0, a);					// bottom right-forward
			f.element[f.index(Nx, Ny, -1, a)] = f(Nx - 1, Ny - 1, 0, a);			// bottom right-back
		}
		//-----------------------------------------------------------------------------------------------------------
		// set left side ghost points
		if (BCtype[0] == 1) {	// Periodic BC
			ix = -1;
			for (iz = 0; iz < Nz; iz++) {
				f(-1, -1, iz, a) = f(Nx - 1, Ny - 1, iz, a);	// set left-forward corner side
				f(-1, Ny, iz, a) = f(Nx - 1, 0, iz, a);			// set left-back corner side
				for (iy = 0; iy < Ny; iy++) {
					f.element[f.index(ix, iy, iz, a)] = f(ix + Nx, iy, iz, a);
				}
			}
		}
		if (BCtype[0] == 2 && t == -1) {	// Velocity BC, inlet
			ix = -1;
			for (iz = 0; iz < Nz; iz++) {
				f(-1, -1, iz, a) = calcfeq(rho, e, Bvel, 0, 0, 0, a, 0);			// set left-forward edge 
				f(-1, Ny, iz, a) = calcfeq(rho, e, Bvel, 0, 0, 0, a, 0);			// set left-back edge
				for (iy = 0; iy < Ny; iy++) {
					f.element[f.index(ix, iy, iz, a)] = calcfeq(rho, e, Bvel, 0, 0, 0, a, 0);
				}
			}

			// copying to left upper edge
			for (iy = 0; iy < Ny; iy++)
				f.element[f.index(-1, iy, Nz, a)] = calcfeq(rho, e, Bvel, Nx - 1, iy, 0, a, 0);
			// copying to left lower edge
			for (iy = 0; iy < Ny; iy++)
				f.element[f.index(-1, iy, -1, a)] = calcfeq(rho, e, Bvel, 0, 0, 0, a, 0);

			f.element[f.index(-1, -1, Nz, a)] = calcfeq(rho, e, Bvel, Nx - 1, Ny - 1, 0, a, 0);		// upper left-forward
			f.element[f.index(-1, Ny, Nz, a)] = calcfeq(rho, e, Bvel, Nx - 1, 0, 0, a, 0);				// upper left-back
			f.element[f.index(-1, -1, -1, a)] = calcfeq(rho, e, Bvel, Nx - 1, Ny - 1, Nz - 1, a, 0);	// bottom left-forward
			f.element[f.index(-1, Ny, -1, a)] = calcfeq(rho, e, Bvel, Nx - 1, 0, Nz - 1, a, 0);		// bottom left-back

		}
		//-----------------------------------------------------------------------------------------------------------
		// set back side ghost points
		if (BCtype[3] == 1) {	// Periodic BC
			iy = Ny;
			for (iz = 0; iz < Nz; iz++) {
				for (ix = 0; ix < Nx; ix++) {
					f.element[f.index(ix, iy, iz, a)] = f(ix, 0, iz, a);
				}
			}
		}
		if (BCtype[3] == 2 && t == -1) {	// Velocity BC, outlet
			iy = Ny;
			for (iz = 0; iz < Nz; iz++) {
				for (ix = 0; ix < Nx; ix++) {
					f.element[f.index(ix, iy, iz, a)] = f.element[f.index(ix, iy - 1, iz, a)];
				}
			}
		}
		// set front side ghost points
		if (BCtype[2] == 1) {	// Periodic BC
			iy = -1;
			for (iz = 0; iz < Nz; iz++) {
				for (ix = 0; ix < Nx; ix++) {
					f.element[f.index(ix, iy, iz, a)] = f(ix, iy + Ny, iz, a);
				}
			}
		}
		if (BCtype[2] == 2 && t == -1) {	// Velocity BC, inlet
			iy = -1;
			for (iz = 0; iz < Nz; iz++) {
				for (ix = 0; ix < Nx; ix++) {
					f.element[f.index(ix, iy, iz, a)] = calcfeq(rho, e, Bvel, 1, 1, 1, a, 0);
				}
			}
		}

	}
}

void updatePBC(direction_density& f) {
	int a = 0;
	int ix = 0;
	int iy = 0;
	int iz = 0;
	//	switch (choice) {
	//	case 1:
	for (a = 0; a < 27; a++) {
		//----------------------------------------------------------------------
		// upper boundary layer (copying from real bottom layer)
		// main surface
		iz = 0;
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				f.element[f.index(ix, iy, iz + Nz, a)] = f(ix, iy, iz, a);
			}
		}
		// edges
		// copying to left upper edge
		for (iy = 0; iy < Ny; iy++)
			f.element[f.index(-1, iy, iz + Nz, a)] = f(Nx - 1, iy, iz, a);

		// copying to right upper edge
		for (iy = 0; iy < Ny; iy++)
			f.element[f.index(Nx, iy, iz + Nz, a)] = f(0, iy, iz, a);

		// copying to forward upper edge
		for (ix = 0; ix < Nx; ix++)
			f.element[f.index(ix, -1, iz + Nz, a)] = f(ix, Ny - 1, iz, a);

		// copying to back upper edge
		for (ix = 0; ix < Nx; ix++)
			f.element[f.index(ix, Ny, iz + Nz, a)] = f(ix, 0, iz, a);

		// corners
		f.element[f.index(-1, -1, Nz, a)] = f(Nx - 1, Ny - 1, 0, a);		// upper left-forward
		f.element[f.index(Nx, -1, Nz, a)] = f(0, Ny - 1, 0, a);			// upper right-forward
		f.element[f.index(Nx, Ny, Nz, a)] = f(0, 0, 0, a);	// upper right-back
		f.element[f.index(-1, Ny, Nz, a)] = f(Nx - 1, 0, 0, a);		// upper left-back
																	//----------------------------------------------------------------------
																	// lower boundary layer (copiyng from real upper layer)
																	// main surface
		iz = Nz - 1;
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				f.element[f.index(ix, iy, -1, a)] = f(ix, iy, iz, a);
			}
		}
		// edges
		// copying to left lower edge
		for (iy = 0; iy < Ny; iy++)
			f.element[f.index(-1, iy, -1, a)] = f(Nx - 1, iy, iz, a);

		// copying to right lower edge
		for (iy = 0; iy < Ny; iy++)
			f.element[f.index(Nx, iy, -1, a)] = f(0, iy, iz, a);

		// copying to forward lower edge
		for (ix = 0; ix < Nx; ix++)
			f.element[f.index(ix, -1, -1, a)] = f(ix, Ny - 1, iz, a);

		// copying to back lower edge
		for (ix = 0; ix < Nx; ix++)
			f.element[f.index(ix, Ny, -1, a)] = f(ix, 0, iz, a);

		// corners
		f.element[f.index(-1, -1, -1, a)] = f(Nx - 1, Ny - 1, Nz - 1, a);		// bottom left-forward
		f.element[f.index(Nx, -1, -1, a)] = f(0, Ny - 1, Nz - 1, a);		// bottom right-forward
		f.element[f.index(Nx, Ny, -1, a)] = f(0, 0, Nz - 1, a);			// bottom right-back
		f.element[f.index(-1, Ny, -1, a)] = f(Nx - 1, 0, Nz - 1, a);		// bottom left-ba
																			//----------------------------------------------------------------------
																			// copy to right side
		ix = Nx;
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				f.element[f.index(ix, iy, iz, a)] = f(0, iy, iz, a);
			}
		}
		// copy to left side
		ix = -1;
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				f.element[f.index(ix, iy, iz, a)] = f(ix + Nx, iy, iz, a);
			}
		}
		// copy to back side
		iy = Ny;
		for (iz = 0; iz < Nz; iz++) {
			for (ix = 0; ix < Nx; ix++) {
				f.element[f.index(ix, iy, iz, a)] = f(ix, 0, iz, a);
			}
		}
		// copy to forward side
		iy = -1;
		for (iz = 0; iz < Nz; iz++) {
			for (ix = 0; ix < Nx; ix++) {
				f.element[f.index(ix, iy, iz, a)] = f(ix, iy + Ny, iz, a);
			}
		}
		for (iz = 0; iz < Nz; iz++) {
			f(-1, -1, iz, a) = f(Nx - 1, Ny - 1, iz, a);	// copy to left-forward corner side
			f(Nx, -1, iz, a) = f(0, Ny - 1, iz, a);			// copy to right-forward corner side
			f(-1, Ny, iz, a) = f(Nx - 1, 0, iz, a);			// copy to left-back corner side
			f(Nx, Ny, iz, a) = f(0, 0, iz, a);				// copy to right-back corner side
		}
	}
}

void updatePBC_solid(Solid_list& solid_list) {
	int ix = 0;
	int iy = 0;
	int iz = 0;
	//----------------------------------------------------------------------
	// upper boundary layer (copying from real bottom layer)
	// main surface
	iz = 0;
	for (iy = 0; iy < Ny; iy++) {
		for (ix = 0; ix < Nx; ix++) {
			solid_list(ix, iy, iz + Nz) = solid_list(ix, iy, iz);
		}
	}
	// edges
	// copying to left upper edge
	for (iy = 0; iy < Ny; iy++)
		solid_list(-1, iy, iz + Nz) = solid_list(Nx - 1, iy, iz);

	// copying to right upper edge
	for (iy = 0; iy < Ny; iy++)
		solid_list(Nx, iy, iz + Nz) = solid_list(0, iy, iz);

	// copying to forward upper edge
	for (ix = 0; ix < Nx; ix++)
		solid_list(ix, -1, iz + Nz) = solid_list(ix, Ny - 1, iz);

	// copying to back upper edge
	for (ix = 0; ix < Nx; ix++)
		solid_list(ix, Ny, iz + Nz) = solid_list(ix, 0, iz);

	// corners
	solid_list(-1, -1, Nz) = solid_list(Nx - 1, Ny - 1, 0);		// upper left-forward
	solid_list(Nx, -1, Nz) = solid_list(0, Ny - 1, 0);			// upper right-forward
	solid_list(Nx, Ny, Nz) = solid_list(0, 0, 0);	// upper right-back
	solid_list(-1, Ny, Nz) = solid_list(Nx - 1, 0, 0);		// upper left-back
															//----------------------------------------------------------------------
															// lower boundary layer (copiyng from real upper layer)
															// main surface
	iz = Nz - 1;
	for (iy = 0; iy < Ny; iy++) {
		for (ix = 0; ix < Nx; ix++) {
			solid_list(ix, iy, -1) = solid_list(ix, iy, iz);
		}
	}
	// edges
	// copying to left lower edge
	for (iy = 0; iy < Ny; iy++)
		solid_list(-1, iy, -1) = solid_list(Nx - 1, iy, iz);

	// copying to right lower edge
	for (iy = 0; iy < Ny; iy++)
		solid_list(Nx, iy, -1) = solid_list(0, iy, iz);

	// copying to forward lower edge
	for (ix = 0; ix < Nx; ix++)
		solid_list(ix, -1, -1) = solid_list(ix, Ny - 1, iz);

	// copying to back lower edge
	for (ix = 0; ix < Nx; ix++)
		solid_list(ix, Ny, -1) = solid_list(ix, 0, iz);

	// corners
	solid_list(-1, -1, -1) = solid_list(Nx - 1, Ny - 1, Nz - 1);		// bottom left-forward
	solid_list(Nx, -1, -1) = solid_list(0, Ny - 1, Nz - 1);		// bottom right-forward
	solid_list(Nx, Ny, -1) = solid_list(0, 0, Nz - 1);			// bottom right-back
	solid_list(-1, Ny, -1) = solid_list(Nx - 1, 0, Nz - 1);		// bottom left-ba
																//----------------------------------------------------------------------
																// copy left side
	ix = Nx;
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			solid_list(ix, iy, iz) = solid_list(0, iy, iz);
		}
	}
	// copy right side
	ix = -1;
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			solid_list(ix, iy, iz) = solid_list(ix + Nx, iy, iz);
		}
	}
	// copy forward side
	iy = Ny;
	for (iz = 0; iz < Nz; iz++) {
		for (ix = 0; ix < Nx; ix++) {
			solid_list(ix, iy, iz) = solid_list(ix, 0, iz);
		}
	}
	// copy back side
	iy = -1;
	for (iz = 0; iz < Nz; iz++) {
		for (ix = 0; ix < Nx; ix++) {
			solid_list(ix, iy, iz) = solid_list(ix, iy + Ny, iz);
		}
	}
	for (iz = 0; iz < Nz; iz++) {
		solid_list(-1, -1, iz) = solid_list(Nx - 1, Ny - 1, iz);	// copy left-forward corner side
		solid_list(Nx, -1, iz) = solid_list(0, Ny - 1, iz);			// copy right-forward corner side
		solid_list(-1, Ny, iz) = solid_list(Nx - 1, 0, iz);			// copy left-back corner side
		solid_list(Nx, Ny, iz) = solid_list(0, 0, iz);				// copy right-back corner side
	}
}

void printstuff(FILE * velfile, FILE * densfile, FILE * parfile, FILE * reyfile, FILE * stressfile, FILE * forcefile, FILE * nhatfile, FILE * sttensfile, FILE * torfile, FILE * erodefile, FILE * eronumbfile, FILE * dmfile, FILE * eroforcefile, FILE * volumefile, FILE * surfacefile, int t, velocity& u, density& rho, vector3Ncubed& tau_stress, vector3Ncubed& F_D, Normalvector& nhat, Stresstensor& stresstensor, vector3Ncubed& torque, density& erodelist, vectorNcubed& F_vdw, Solid_list& solid_list, density& masschange) {
	int ix = 0;
	int iy = 0;
	int iz = 0;
	double F[3] = { 0 };
	double F_vdwsq = 0;
	double F_stress = 0;
	int vol = 0;
	int sur = 0;
//	if (t == 0) {
//		cout << "\n tau = " << tau << "\n";
//		cout << "\n Lz = " << Lz << "\n";
//		cout << "\n gg = " << gg << "\n";
//		cout << "\n F = " << F << "\n";
//		fprintf(parfile, "%i %i %i %e %e %e %e %e %e %e", Nx, Ny, Nz, Re, tau, Delta_T, sphere_radius, masspernode, kappa_er, VDW_0);
//	}
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				fprintf(velfile, "%e %e %e ", u(ix, iy, iz, 0), u(ix, iy, iz, 1), u(ix, iy, iz, 2));
				//fprintf(densfile, "%e ", rho(ix, iy, iz));
				fprintf(stressfile, "%e %e %e ", tau_stress(ix, iy, iz, 0), tau_stress(ix, iy, iz, 1), tau_stress(ix, iy, iz, 2));
				fprintf(forcefile, "%e ", F_D(ix, iy, iz, 0));
				fprintf(nhatfile, "%e %e %e ", nhat(ix, iy, iz, 0), nhat(ix, iy, iz, 1), nhat(ix, iy, iz, 2));
				fprintf(torfile, "%e %e %e ", torque(ix, iy, iz, 0), torque(ix, iy, iz, 1), torque(ix, iy, iz, 2));
				//fprintf(erodefile, "%i ", erodelist(ix, iy, iz));
				//fprintf(dmfile, "%e ", masschange(ix, iy, iz));
				if (solid_list(ix, iy, iz) == 0) {
					F_stress = sqrt(pow(tau_stress(ix, iy, iz, 0), 2) + pow(tau_stress(ix, iy, iz, 1), 2) + pow(tau_stress(ix, iy, iz, 2), 2));
					fprintf(eroforcefile, "%e ", F_stress / F_vdw(ix, iy, iz));
					fprintf(eronumbfile, "%e ", kappa_er*dt*(F_stress - F_vdw(ix, iy, iz)) / masspernode);
				}
				else if (solid_list(ix, iy, iz) == -1 || solid_list(ix, iy, iz) == 1) {
					fprintf(eroforcefile, "%e ", 0.);
					fprintf(eronumbfile, "%e ", 0.);
				}
				//F[0] += F_D(ix, iy, iz, 0);
				//F[1] += F_D(ix, iy, iz, 1);
				//F[2] += F_D(ix, iy, iz, 2);
			}
		}
	}
	/*	if (t == 3000) {
	for (int a = 0; a < 3; a++) {
	for (int b = 0; b < 3; b++) {

	for (iz = 0; iz < Nz; iz++) {
	for (iy = 0; iy < Ny; iy++) {
	for (ix = 0; ix < Nx; ix++) {
	fprintf(sttensfile, "%e ", stresstensor(ix, iy, iz, a, b));
	}
	}
	}
	fprintf(sttensfile, "\n");
	}
	}
	}*/
	fprintf(forcefile, "\n");
	//fprintf(forcefile, "%e %e %e\n", F[0], F[1], F[2]);
	fprintf(velfile, "\n");
	fprintf(densfile, "\n");
	fprintf(stressfile, "\n");
	fprintf(nhatfile, "\n");
	fprintf(torfile, "\n");
	fprintf(erodefile, "\n");
	fprintf(eronumbfile, "\n");
	fprintf(dmfile, "\n");
	fprintf(eroforcefile, "\n");
}

vector<int> readBC(string x0BC, string xNBC, string y0BC, string yNBC, string z0BC, string zNBC) {
	int ix = 0;
	int iy = 0;
	int iz = 0;
	int a = 0;
	int i = 0;
	vector<string> BC = { "periodicBC", "velocityBC" };	// if an additional boundary is implemented, loop the for loop below iend + 1 times.
	vector<int> BCtype;
	BCtype.resize(6);
	fill(BCtype.begin(), BCtype.end(), 0);

	for (i = 0; i < 2; i++) {
		if (x0BC.compare(BC[i]) == 0)
			BCtype[0] = 1 + i;
		if (xNBC.compare(BC[i]) == 0)
			BCtype[1] = 1 + i;
		if (y0BC.compare(BC[i]) == 0)
			BCtype[2] = 1 + i;
		if (yNBC.compare(BC[i]) == 0)
			BCtype[3] = 1 + i;
		if (z0BC.compare(BC[i]) == 0)
			BCtype[4] = 1 + i;
		if (zNBC.compare(BC[i]) == 0)
			BCtype[5] = 1 + i;
	}
	for (i = 0; i < 6; i++) {
		if (BCtype[i] == 0)
			cout << "Error in BC. Edge " << i << " has no BC. Assign a BC to this edge.";
	}

	return BCtype;



}

double find_umax(velocity& u) {
	int iz = 0;
	int iy = 0;
	int ix = 0;
	double umax[3] = { 0 };
	double umaxtemp[3] = { 0 };
	double ans[4];
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				umaxtemp[0] = u(ix, iy, iz, 0);
				umaxtemp[1] = u(ix, iy, iz, 1);
				umaxtemp[2] = u(ix, iy, iz, 2);
				if (abs(umaxtemp[0]) > umax[0]) {
					umax[0] = umaxtemp[0];

				}
				if (abs(umaxtemp[1]) > umax[1])
					umax[1] = umaxtemp[1];
				if (abs(umaxtemp[2]) > umax[2])
					umax[2] = umaxtemp[2];
			}
		}
	}
	//	ans[0] = sqrt(umax[0] * umax[0] + umax[1] * umax[1] + umax[2] * umax[2]);
	//	ans[1] = 
	return sqrt(umax[0] * umax[0] + umax[1] * umax[1] + umax[2] * umax[2]);

}
double find_uav(velocity& u) {
	int ix = 0;
	int iy = 0;
	int iz = 0;
	double uxsum = 0;
	double uysum = 0;
	double uzsum = 0;
	double uav;

	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				uxsum += u(ix, iy, iz, 0);
				uysum += u(ix, iy, iz, 1);
				uzsum += u(ix, iy, iz, 2);
			}
		}
	}
	uxsum = uxsum / ((double)Nx*(double)Ny*(double)Nz);
	uysum = uysum / ((double)Nx*(double)Ny*(double)Nz);
	uzsum = uzsum / ((double)Nx*(double)Ny*(double)Nz);
	uav = sqrt(uxsum*uxsum + uysum*uysum + uzsum*uzsum);
	return uav;
}
