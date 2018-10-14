#include "LBM_fluidsteps_func.h"
#include "LBM_force_erosion_func.h"


void IC(momentum_direction& e, direction_density& f, direction_density& ftemp, Solid_list& solid_list, vectorNcubed& F_vdw) {
	// variables
	// -----------------------------
	int a = 0;
	int ix = 0;		// x 
	int iy = 0;		// y
	int iz = 0;		// z
	int cellist[27] = { 3, 2, 3, 2, 1, 2, 3, 2, 3  ,  2, 1, 2, 1, 0, 1, 2, 1, 2  ,  3, 2, 3, 2, 1, 2, 3, 2, 3 };
	double weights[4];
	double w = 0;
	weights[0] = 0.296296296296296;//8. / 27.;	
	weights[1] = 0.074074074074074;//2. / 27.;
	weights[2] = 0.018518518518519;// 1. / 54.;
	weights[3] = 0.004629629629630;//1. / 216.;
								   //------------------------------
								   // function body
								   //		cout << "\n IC: ------------------------------------" << "\n";
	for (a = 0; a < 27; a++) {
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					if (solid_list(ix, iy, iz) != 1) {
						w = weights[cellist[a]];
						f(ix, iy, iz, a) = w;
						ftemp(ix, iy, iz, a) = w;
					}
					else {
						f(ix, iy, iz, a) = 0.;
						ftemp(ix, iy, iz, a) = 0.;
					}
					if (a == 0) {
						computeVDWforce(ix, iy, iz, e, solid_list, F_vdw);
					}

				}
			}
		}
	}

}


void macrovariables(velocity& u, density& rho, Solid_list& solid_list, direction_density& f, momentum_direction& e) {
	// variables
	int ix = 0;
	int iy = 0;
	int iz = 0;
	int a = 0;		// momentum direction
					// function body

					//	cout << "\n Macro: ------------------------------------ \n";
	for (a = 0; a < 27; a++) {   //Could maybe swap 27 for 3^3, where 3 is dimension.
								 //		cout << "\n a = " << a << "\n";
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					//					cout << " " << f(ix, iy, iz, a) << " ";
					if (a == 0) {
						u(ix, iy, iz, 0) = 0.;
						u(ix, iy, iz, 1) = 0.;
						u(ix, iy, iz, 2) = 0.;
						rho(ix, iy, iz) = 0.;
					}

					if (solid_list(ix, iy, iz) == -1) {
						rho(ix, iy, iz) += f(ix, iy, iz, a);
						u(ix, iy, iz, 0) += e(a, 0) * f(ix, iy, iz, a);
						u(ix, iy, iz, 1) += e(a, 1) * f(ix, iy, iz, a);
						u(ix, iy, iz, 2) += e(a, 2) * f(ix, iy, iz, a);
					}
					//					cout << "x: " << u(ix,iy,iz,0) << ", y: " << u(ix,iy,iz,1) << ", z: " << u(ix,iy,iz,2) << "\n";
					//					cout << " " << e(a,0) << " " << e(a,1) << " " << e(a,2);
					//					if (ix == 1 && iy == 1 && iz == 0)
					//						cout << " " << f(ix, iy, iz, a) << " ";
					//					if (u(ix, iy, iz, 0) > 5000 || u(ix, iy, iz, 1) > 5000 || u(ix, iy, iz, 2) > 5000)
					//						cout << "Warning: Very high velocity!! \n";
				}
			}
			//			cout << "\n";
		}
	}
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				if (solid_list(ix, iy, iz) == -1) {
					//					cout << " " << u(ix, iy, iz, 0) << " " << u(ix, iy, iz, 1) << " " << u(ix, iy, iz, 2) << " " << rho(ix, iy, iz) << "\n";
					if (rho(ix, iy, iz) == 0.) {
						cout << "Error: rho = 0\n";
						break;
					}
					u(ix, iy, iz, 0) /= rho(ix, iy, iz);
					u(ix, iy, iz, 1) /= rho(ix, iy, iz);
					u(ix, iy, iz, 2) /= rho(ix, iy, iz);
				}
			}
		}
	}


};



void edf(Solid_list& solid_list, velocity& u, density& rho, EDF& feq, momentum_direction& e, int forcedirection) {
	// variables
	double f1 = 3.;
	double f2 = 9.*0.5;
	double f3 = 3.*0.5;
	double w = 0.;
	double usq = 0.;
	double ueq[3];
	double csq = 0.;
	double dotprod = 0.;
	int ix = 0;
	int iy = 0;
	int iz = 0;
	int a = 0;	// Direction counter. Used for feq and w. 
				// function body
				// cout << "\n EDF: -------------------------------------- \n";
	for (a = 0; a < 27; a++) {
		//		cout << "\n a = " << a << ": \n \n";
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {

					if (solid_list(ix, iy, iz) == -1) {
						/*						ueq[0] = u(ix, iy, iz, 0) + (tau*gg) / rho(ix, iy, iz);	//Add subroutine for this;
						ueq[1] = u(ix, iy, iz, 1);
						ueq[2] = u(ix, iy, iz, 2);
						w = weights[cellist[a]];
						//							cout << "\n " << ux << " " << uy << " " << uz << "\n";
						usq = ueq[0] * ueq[0] + ueq[1] * ueq[1] + ueq[2] * ueq[2];
						csq = c*c;
						dotprod = edfvecdot(e, ueq, a);		// Testa så att beräkna feq i den nya funktionen funkar.
						feq(ix, iy, iz, a) = rho(ix, iy, iz) * w * (1. + f1*dotprod / (csq)+f2*dotprod*dotprod / (csq*csq) - f3*usq / csq);*/ //c, c^2, c^2
						ueq[0] = u(ix, iy, iz, 0);
						ueq[1] = u(ix, iy, iz, 1);
						ueq[2] = u(ix, iy, iz, 2);
						feq(ix, iy, iz, a) = calcfeq(rho, e, ueq, ix, iy, iz, a, forcedirection);
					}
					//					if (ix == 1 && iy == 1 && iz == 0)
					//						cout << " " << feq(ix, iy, iz, a) << " ";
				}
			}
			//			cout << "\n";
		}
		//			cout << "\n";
	}
}
double edfvecdot(momentum_direction& e, double ueq[3], int a) {
	double answer = 0;
	answer = (double)e(a, 0) * ueq[0] + (double)e(a, 1) * ueq[1] + (double)e(a, 2) * ueq[2];
	return answer;
};
double calcfeq(density& rho, momentum_direction& e, double ueq[3], int ix, int iy, int iz, int a, int forcedirection) {
	double weights[4] = { 0 };
	weights[0] = 0.296296296296296;// 8. / 27.;		
	weights[1] = 0.074074074074074;// 2. / 27.;
	weights[2] = 0.018518518518519;// 1. / 54.;
	weights[3] = 0.004629629629630;// 1. / 216.;
	double feqval;
	double f1 = 3.;
	double f2 = 9.*0.5;
	double f3 = 3.*0.5;
	double w = 0.;
	double usq = 0.;
	double csq = 0.;
	double dotprod = 0.;
	ueq = calcforce(ueq, rho, ix, iy, iz, forcedirection);
	w = weights[cellist[a]];
	//							cout << "\n " << ux << " " << uy << " " << uz << "\n";
	usq = ueq[0] * ueq[0] + ueq[1] * ueq[1] + ueq[2] * ueq[2];
	csq = c*c;
	dotprod = (double)e(a, 0) * ueq[0] + (double)e(a, 1) * ueq[1] + (double)e(a, 2) * ueq[2];
	feqval = rho(ix, iy, iz) * w * (1. + f1*dotprod / (csq)+f2*dotprod*dotprod / (csq*csq) - f3*usq / csq);	// c, c^2, c^2.
	return feqval;
}
double * calcforce(double ueq[3], density& rho, int ix, int iy, int iz, int forcedirection) { //forcedirection == "0" betyder ingen kraft verkar på systemet.
	if (forcedirection == 1)
		ueq[0] += tau*gg / rho(ix, iy, iz); //F*cos(((double)iz / Lz)*pi);	//tau*gg / rho(ix, iy, iz);
	if (forcedirection == 2)
		ueq[1] += F*cos(((double)iz / Lz)*pi);  //tau*gg / rho(ix, iy, iz);
	if (forcedirection == 3)
		ueq[2] += F*cos(((double)iz / Lz)*pi);  //tau*gg / rho(ix, iy, iz);
	return ueq;
}


void stream(Solid_list& solid_list, direction_density& f, direction_density& ftemp, momentum_direction& e) {
	// variable defs
	int ix = 0;			// x
	int iy = 0;			// y
	int iz = 0;			// z
	int a = 0;			// direction counter
	int ashift = 0;
	// function body
	/*	cout << "\n Stream: ------------------------------------ \n";
	for (a = 0; a < 27; a++) {
	cout << "\n a = " << a << "\n";
	for (iz = -1; iz < Nz + 1; iz++) {
	//			cout << "\n";
	for (iy = -1; iy < Ny + 1; iy++) {
	for (ix = -1; ix < Nx + 1; ix++) {
	if (ix == 1 && iy == 0 && iz == 0)
	cout << " " << f(ix, iy, iz, a) << " ";
	}
	}
	}
	}*/

	for (a = 0; a < 27; a++) {
		ashift = 2 * (13 - a);
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					switch (solid_list(ix, iy, iz)) {
					case 1:	// interior solid node = Do nothing
						break;
					case 0:	// surface solid node = bounceback BC
						ftemp(ix, iy, iz, a) = f(ix + e(26 - a, 0), iy + e(26 - a, 1), iz + e(26 - a, 2), a);
						break;
					case -1: // fluid node = regular stream
						ftemp(ix, iy, iz, a) = f(ix + e(26 - a, 0), iy + e(26 - a, 1), iz + e(26 - a, 2), a);	//the e:s shifts ix,iy,iz to the node streaming into the node we're looking at. 
						break;
					}
				}
			}
		}
	}

	/*	cout << "---------------------------------------------------------------";
	for (a = 0; a < 27; a++) {
	cout << "\n a = " << a << "\n";
	for (iz = -1; iz < Nz + 1; iz++) {
	//			cout << "\n";
	for (iy = -1; iy < Ny + 1; iy++) {
	for (ix = -1; ix < Nx + 1; ix++) {
	if (ix == 1 && iy == 0 && iz == 0)
	cout << " " << ftemp(ix, iy, iz, a) << " ";
	}
	}
	}
	}*/
}
void bouncebackBC(direction_density& ftemp, direction_density& f, momentum_direction& e, int ix, int iy, int iz, int a, int ashift) {

	/*	int ixtemp = 0;
	int iytemp = 0;
	int iztemp = 0;
	ixtemp = ix + e(26 - a, 0);
	iytemp = iy + e(26 - a, 1);
	iztemp = iz + e(26 - a, 2);
	if (ixtemp < 0)
	ixtemp += Nx;
	if (ixtemp >(Nx - 1))
	ixtemp -= Nx;
	if (iytemp < 0)
	iytemp += Ny;
	if (iytemp >(Ny - 1))
	iytemp -= Ny;
	if (iztemp < 0)
	iztemp += Nz;
	if (iztemp >(Nz - 1))
	iztemp -= Nz;*/
	//	cout << ftemp(ixtemp, iytemp, iztemp, a + ashift) << " " << f(ixtemp, iytemp, iztemp, a) << "\n"; 
	ftemp(ix, iy, iz, a + ashift) = f(ix, iy, iz, a);
	//	cout << ftemp(ixtemp, iytemp, iztemp, a + ashift) << " " << f(ixtemp, iytemp, iztemp, a);
}


void collision(Solid_list& solid_list, direction_density& f, direction_density& ftemp, EDF& feq, momentum_direction& e) {
	int ix = 0;
	int iy = 0;
	int iz = 0;

	int ixtemp;
	int iytemp;
	int iztemp;
	int a = 0;
	int ashift = 0;
	//	cout << "\n Collision: ------------------------------------\n \n";
	/*	for (a = 0; a < 27; a++) {
	cout << "\n a = " << a << "\n";
	for (iz = 0; iz < Nz; iz++) {
	cout << "\n";
	for (iy = 0; iy < Ny; iy++) {
	for (ix = 0; ix < Nx; ix++) {
	if (ix == 1 && iy == 1 && iz == 0)
	cout << " " << ftemp(ix, iy, iz, a) << " ";
	}
	}
	}
	}*/

	for (a = 0; a < 27; a++) {
		ashift = 2 * (13 - a);
		//		cout << "\n a = " << a << "\n";
		for (iz = 0; iz < Nz; iz++) {
			//				cout << "\n";
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					switch (solid_list(ix, iy, iz)) {
					case 1:	// Interior node
						break;
					case 0:	// Surface node
						f(ix, iy, iz, a + ashift) = ftemp(ix, iy, iz, a);
						break;
					case -1: //	Fluid node
						f(ix, iy, iz, a) = ftemp(ix, iy, iz, a) - (ftemp(ix, iy, iz, a) - feq(ix, iy, iz, a)) / tau;
						break;
					}
				}
			}
			//			cout << "\n";
		}
	}

	/*	for (a = 0; a < 27; a++) {
	cout << "\n a = " << a << "\n";
	for (iz = 0; iz < Nz; iz++) {
	cout << "\n";
	for (iy = 0; iy < Ny; iy++) {
	for (ix = 0; ix < Nx; ix++) {
	if (ix == 1 && iy == 1 && iz == 0)
	cout << " " << f(ix, iy, iz, a) << " ";
	}
	}
	}
	}*/
};