#include "LBM_force_erofunc.h"

void computestress(momentum_direction& e, direction_density& ftemp, direction_density& f, EDF& feq, Solid_list& solid_list, Stresstensor& stresstensor, Normalvector& nhat, Wall_force& tau_stress, vector3Ncubed& F_sum) {
	int ix = 0;
	int iy = 0;
	int iz = 0;
	int a = 0;
	int i = 0;	//e component index. Row index. i:th component of the stress.
	int j = 0;	//e component index. Column index. i:th component of the stress acting on the j:th surface.
	int ixtemp = 0; //temp nhat x
	int iytemp = 0; //temp nhat y
	int iztemp = 0; //temp nhat z
	int normcount = 0;
	double normvec[3] = { 0. };
	double normfactvec[3] = { sq3inv, sq2inv, 1. };
	double normfactor = 0;
	int Nf = 0; //number of nnn fluid points

				// Calculating stress tensor and normal vector
				/*	cout << "upper: \n";
				for (a = 0; a < 27; a++) {
				cout << " " << f(1, 5, 1, a) << " ";
				}
				cout << "\n lower: \n";
				for (a = 0; a < 27; a++) {
				cout << " " << f(1, 0, 1, a) << " ";
				}*/
	stresstensor.clear();
	for (a = 0; a < 27; a++) {
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {

					if (solid_list(ix, iy, iz) == -1) { // Fluid node
						for (i = 0; i < 3; i++) {
							for (j = 0; j < 3; j++) {
								stresstensor(ix, iy, iz, i, j) += (1. - (1. / (2.*tau)))*(ftemp(ix, iy, iz, a) - feq(ix, iy, iz, a))*e(a, i)*e(a, j);
							}
						}
					}
					if (solid_list(ix, iy, iz) == 0) { // Surface solid node
						if (solid_list(ix + e(a, 0), iy + e(a, 1), iz + e(a, 2)) != -1) {
							nhat(ix, iy, iz, 0) += e(26 - a, 0);
							nhat(ix, iy, iz, 1) += e(26 - a, 1);
							nhat(ix, iy, iz, 2) += e(26 - a, 2);
						}
					}

				}
			}
		}
	}
	// Normalizing normal vector. We only care about direction and that abs(nhat) = 1. 1/sqrt(3)*(1,1,1) = 1/sqrt(12)*(2,2,2).
	// So we set nhats elements = 1, 0, -1 so we can use them as indeces too. nhat essentially becomes a normalized e vector (mom_dir).
	// This works for our appprox surface. For a real sphere ex. then this is wrong. Then we need nhat = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta)).
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				if (solid_list(ix, iy, iz) == 0) { //Surface solid point
					normcount = 0;
					for (i = 0; i < 3; i++) {
						if (nhat(ix, iy, iz, i) > 0)
							nhat(ix, iy, iz, i) = 1;
						if (nhat(ix, iy, iz, i) < 0)
							nhat(ix, iy, iz, i) = -1;
						if (nhat(ix, iy, iz, i) == 0)
							normcount++;
					}
					if (normcount == 3)
						cout << "\n Error in computestress! normal vector = (0, 0, 0).\n";
					normfactor = normvec[normcount];
					/*					ixtemp = nhat(ix, iy, iz, 0);
					iytemp = nhat(ix, iy, iz, 1);
					iztemp = nhat(ix, iy, iz, 2);
					normfactor = sqrt(ixtemp*ixtemp + iytemp*iytemp + iztemp*iztemp);
					nhat(ix, iy, iz, 0) /= normfactor;
					nhat(ix, iy, iz, 1) /= normfactor;
					nhat(ix, iy, iz, 2) /= normfactor;*/
				}

			}
		}
	}
	// Computing wall shear stress and drag force (actually the total force, but for now its just drag).
	//F_D.clear();
	tau_stress.clear();
	//	for (a = 0; a < 27; a++) {
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {

				if (solid_list(ix, iy, iz) == 0) { // solid surface point
												   /*						if (solid_list(ix + e(a, 0), iy + e(a, 1), iz + e(a, 2)) == -1) { // fluid point
												   F_D(ix, iy, iz, 0) += e(26 - a, 0)*(-f(ix + e(a, 0), iy + e(a, 1), iz + e(a, 2), 26 - a) + f(ix, iy, iz, a));
												   F_D(ix, iy, iz, 1) += e(26 - a, 1)*(-f(ix + e(a, 0), iy + e(a, 1), iz + e(a, 2), 26 - a) + f(ix, iy, iz, a));
												   F_D(ix, iy, iz, 2) += e(26 - a, 2)*(-f(ix + e(a, 0), iy + e(a, 1), iz + e(a, 2), 26 - a) + f(ix, iy, iz, a));
												   }*/ //nhat = (e1 + e2 + ...)*normfactor. just sum tau_stress over all fluid points and multiply with normfactor
												   //normvec[0] = 0; normvec[1] = 0; normvec[2] = 0;
												   //for (a = 0; a < 27; a++) {
												   //if (solid_list(ix + e(a, 0), iy + e(a, 1), iz + e(a, 2)) == -1) { // fluid point
												   //normvec[0] += e(a, 0); normvec[1] += e(a, 1); normvec[2] += e(a, 2);
					for (i = 0; i < 3; i++) {
						for (j = 0; j < 3; j++) {
							//F_D(ix, iy, iz, i) += stresstensor(ix + nhat(ix, iy, iz, 0), iy + nhat(ix, iy, iz, 1), iz + nhat(ix, iy, iz, 2), i, j); //* normvec[j]; //sum the force contribution from all surface points. sigma_ij * ehat
							tau_stress(ix, iy, iz, i) += -nhat(ix, iy, iz, j)*stresstensor(ix + nhat(ix, iy, iz, 0), iy + nhat(ix, iy, iz, 1), iz + nhat(ix, iy, iz, 2), i, j); //Should be normfactor * surface area exposed to the fluid. But these 2 cancel out, so no contribution from them.
							F_sum(ix, iy, iz, i) += tau_stress(ix, iy, iz, i);
							//tau_stress(ix, iy, iz, i) += -e(a,j)*stresstensor(ix + e(a, 0), iy + e(a, 1), iz + e(a, 2), i, j);
						}
					}
					//}
					//}
					//normfactor = sqrt(normvec[0] * normvec[0] + normvec[1] * normvec[1] + normvec[2] * normvec[2]);
					//for (i = 0; i < 3; i++)
					//tau_stress(ix, iy, iz, i) = tau_stress(ix, iy, iz, i)*normfactor;
				}



			}
		}
	}
	//	}



}

void computetorque(Solid_list& solid_list, Wall_force& tau_stress, vector3Ncubed& torque) {
	double r[3];
	int ix = 0;
	int iy = 0;
	int iz = 0;
	dvec result;
	result.resize(3);
	torque.clear();
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				if (solid_list(ix, iy, iz) == 0) { // surface node
					r[0] = double(ix) - solid_list.center[0];
					r[1] = double(iy) - solid_list.center[1];
					r[2] = double(iz) - solid_list.center[2];
					result = crossproduct(tau_stress, r, ix, iy, iz);
					torque(ix, iy, iz, 0) = result[0];
					torque(ix, iy, iz, 1) = result[1];
					torque(ix, iy, iz, 2) = result[2];
				}
			}
		}
	}


}

// Function should be after collosion.
void erosion(Solid_list& solid_list, momentum_direction& e, vector3Ncubed& F_sum, density& rho, direction_density& f, int forcedirection, FILE * solfile, density& erodelist, Normalvector& nhat) {
	int ix = 0;
	int iy = 0;
	int iz = 0;
	int ixshift = 0;
	int iyshift = 0;
	int izshift = 0;
	int a = 0;
	double WDWsqsum = 0.; //Wan-der-waals force squared sum.
	double FFsq = 0.; //Fluid Force squared.
	double ueq[3] = { 0 };
	double weights[4] = { 0 };
	weights[0] = 0.296296296296296;// 8. / 27.;	
	weights[1] = 0.074074074074074;// 2. / 27.;
	weights[2] = 0.018518518518519;// 1. / 54.;
	weights[3] = 0.004629629629630;// 1. / 216.;
								   //double Fmaxsq = F_sum.maxsq();


	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				erodelist(ix, iy, iz) = 1.;
				if (solid_list(ix, iy, iz) == 0) { //if we're standing on a surface node
					WDWsqsum = 0.;
					for (a = 0; a < 27; a++) {
						ixshift = ix + e(a, 0);
						iyshift = iy + e(a, 1);
						izshift = iz + e(a, 2);
						if (solid_list(ixshift, iyshift, izshift) == 1 || solid_list(ixshift, iyshift, izshift) == 0) // add wdwforce from all solid nodes.
							WDWsqsum = WDWsqsum + WDWforce;
					}
					FFsq = (pow(F_sum(ix, iy, iz, 0) / (double)updatefreq, 2) + pow(F_sum(ix, iy, iz, 1) / (double)updatefreq, 2) + pow(F_sum(ix, iy, iz, 2) / (double)updatefreq, 2));
					if (FFsq > WDWsqsum) { //if the fluid force is greater than the WDW force from all solid nodes.
						erodelist(ix, iy, iz) = 1;//-kappa_er*sqrt(FFsq - WDWsqsum); //erode point
					}
					else
						erodelist(ix, iy, iz) = 0; //don't erode point
				}
			}
		}
	}
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				if (erodelist(ix, iy, iz) == 1 && solid_list(ix, iy, iz) == 0) {



					solid_list(ix, iy, iz) = -1; //surface node becomes fluid node.
					rho(ix, iy, iz) = rho(ix + nhat(ix, iy, iz, 0), iy + nhat(ix, iy, iz, 1), iz + nhat(ix, iy, iz, 2)); //fluid node is initialized with same density as the interface node.
																														 // check that rho != 0 at new points.
					for (a = 0; a < 27; a++) { //Initialize new fluid point with same f as interface node. Also, check all nearby nodes. If it's a interior solid node, it becomes a surface.
						f(ix, iy, iz, a) = f(ix + nhat(ix, iy, iz, 0), iy + nhat(ix, iy, iz, 1), iz + nhat(ix, iy, iz, 2), a);
						ixshift = ix + e(a, 0);
						iyshift = iy + e(a, 1);
						izshift = iz + e(a, 2);
						if (solid_list(ixshift, iyshift, izshift) == 1) {
							solid_list(ixshift, iyshift, izshift) = 0;
							f(ixshift, iyshift, izshift, a) = weights[cellist[a]];
							rho(ixshift, iyshift, izshift) = 1.;
						}
					}

				}
			}
		}
	}
	solid_list.printsolid_list(solfile);
	erodelist.clear();
	F_sum.clear();
}
