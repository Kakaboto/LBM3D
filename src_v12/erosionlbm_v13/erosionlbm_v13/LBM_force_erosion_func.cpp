#include "LBM_force_erosion_func.h"



void computestress(momentum_direction& e, direction_density& ftemp, direction_density& f, EDF& feq, Solid_list& solid_list, Stresstensor& stresstensor, Normalvector& nhat, vector3Ncubed& tau_stress, vector3Ncubed& F_sum, density& masschange, vectorNcubed& F_vdw, int i_er, int i_Fvdw, density& rho) {
	int ix = 0;
	int iy = 0;
	int iz = 0;
	int a = 0;
	int i = 0;	//e component index. Row index. i:th component of the stress.
	int j = 0;	//e component index. Column index. i:th component of the stress acting on the j:th surface.
	int ixshift = 0; //temp nhat x
	int iyshift = 0; //temp nhat y
	int izshift = 0; //temp nhat z
	int normcount = 0;
	double normvec[3] = { 0. };
	double normfactvec[3] = { sq3inv, sq2inv, 1. };
	double normfactor = 0;
	int Nf = 0; //number of nnn fluid points
	double WDWsum = 0.;
	double FF = 0;
	double shearforce[3] = { 0 };

	// Calculating stress tensor, normal vectors, forces and erosion.
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
					if (normcount == 3) { // This corresponds to an isolated point. In theory this would just float away with the fluid so we erode away that point. Might be thin ice here.
						cout << "\n Isolated point! normal vector = (0, 0, 0).\n";
					}
					normfactor = normfactvec[normcount];
					/*					ixtemp = nhat(ix, iy, iz, 0);
					iytemp = nhat(ix, iy, iz, 1);
					iztemp = nhat(ix, iy, iz, 2);
					normfactor = sqrt(ixtemp*ixtemp + iytemp*iytemp + iztemp*iztemp);
					nhat(ix, iy, iz, 0) /= normfactor;
					nhat(ix, iy, iz, 1) /= normfactor;
					nhat(ix, iy, iz, 2) /= normfactor;*/
					computeVDWforce(ix, iy, iz, e, solid_list, F_vdw);
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
					for (i = 0; i < 3; i++) {
						for (j = 0; j < 3; j++) {
							//shearforce[i] += -nhat(ix, iy, iz, j)*stresstensor(ix + nhat(ix, iy, iz, 0), iy + nhat(ix, iy, iz, 1), iz + nhat(ix, iy, iz, 2), i, j); //Should be normfactor * surface area exposed to the fluid. But these 2 cancel out, so no contribution from them.
							tau_stress(ix,iy,iz,i) += -nhat(ix, iy, iz, j)*stresstensor(ix + nhat(ix, iy, iz, 0), iy + nhat(ix, iy, iz, 1), iz + nhat(ix, iy, iz, 2), i, j); //Should be normfactor * surface area exposed to the fluid. But these 2 cancel out, so no contribution from them.
						}
					}

					//tau_stress(ix, iy, iz) = sqrt(shearforce[0] * shearforce[0] + shearforce[1] * shearforce[1] + shearforce[2] * shearforce[2] - nhat(ix, iy, iz, 0)*shearforce[0] * nhat(ix, iy, iz, 0)*shearforce[0] - nhat(ix, iy, iz, 1)*shearforce[1] * nhat(ix, iy, iz, 1)*shearforce[1] - nhat(ix, iy, iz, 2)*shearforce[2] * nhat(ix, iy, iz, 2)*shearforce[2]);
					//FF = sqrt(pow(tau_stress(ix, iy, iz, 0), 2) + pow(tau_stress(ix, iy, iz, 1), 2) + pow(tau_stress(ix, iy, iz, 2), 2));
					//if (tau_stress(ix,iy,iz) > F_vdw(ix, iy, iz)) { //if the fluid force is greater than the WDW force from all solid nodes.
					//	masschange(ix, iy, iz) += -kappa_er*dt*(tau_stress(ix,iy,iz) - F_vdw(ix,iy,iz)); 

					//random change s� jag kan pusha stuff
					//random change igen....
					FF = sqrt(pow(tau_stress(ix, iy, iz, 0), 2) + pow(tau_stress(ix, iy, iz, 1), 2) + pow(tau_stress(ix, iy, iz, 2), 2));
					if (FF > F_vdw(ix, iy, iz)) { //if the fluid force is greater than the WDW force from all solid nodes.
						masschange(ix, iy, iz) += -kappa_er*dt*(FF - F_vdw(ix,iy,iz)); 
					}
					//else
					//	masschange(ix, iy, iz) += 0.; //don't erode point
				}



			}
		}
	}


}
void computeVDWforce(int ix, int iy, int iz, momentum_direction& e, Solid_list& solid_list, vectorNcubed& F_vdw) {
	double WDWsum = 0.;
	int ixshift = 0;
	int iyshift = 0;
	int izshift = 0;
	for (int a = 0; a < 27; a++) {
		ixshift = ix + e(a, 0);
		iyshift = iy + e(a, 1);
		izshift = iz + e(a, 2);
		if (solid_list(ixshift, iyshift, izshift) == 1 || solid_list(ixshift, iyshift, izshift) == 0) // add wdwforce from all solid nodes.
			WDWsum = WDWsum + VDWforce(e(a, 0), e(a, 1), e(a, 2));
	}
	F_vdw(ix, iy, iz) = WDWsum;
}
double VDWforce(int x, int y, int z) {
	// VDW force scaling with 1/r^2.
	if (x != 0 && y != 0 && z != 0)
		return 0.33333333333333*VDW_0; // r = sqrt(3). 1/3 ~= 0.33333... 
	if ((x != 0 && y != 0) || (x != 0 && z != 0) || (y != 0 && z != 0))
		return 0.5*VDW_0; // r = sqrt(2) => 1/2
	if (x != 0 || y != 0 || z != 0)
		return VDW_0; // r = 1.
					  //if none of the ifs trigger, x = y = z = 0.
	return 0; // the solid node itself.
}

void computetorque(Solid_list& solid_list, vector3Ncubed& tau_stress, vector3Ncubed& torque) {
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
dvec crossproduct(vector3Ncubed& tau_stress, double r[3], int ix, int iy, int iz) {
	// assumed vector of dimension 3.
	dvec result;
	result.resize(3);

	result[0] = r[1] * tau_stress(ix, iy, iz, 2) - r[2] * tau_stress(ix, iy, iz, 1);
	result[1] = -r[0] * tau_stress(ix, iy, iz, 2) + r[2] * tau_stress(ix, iy, iz, 0);
	result[2] = r[0] * tau_stress(ix, iy, iz, 1) - r[1] * tau_stress(ix, iy, iz, 0);

	return result;
}

// Function should be after collosion.
void erosion(Solid_list& solid_list, momentum_direction& e, vector3Ncubed& F_sum, density& rho, velocity& u, direction_density& ftemp, EDF& feq, int forcedirection, FILE * solfile, density& masschange, Normalvector& nhat, vectorNcubed& ero_reso_check, vectorNcubed& F_vdw, FILE * errorfile, int* ecp) {
	int ix = 0;
	int iy = 0;
	int iz = 0;
	int ixshift = 0;
	int iyshift = 0;
	int izshift = 0;
	int a = 0;
	double FFsq = 0.; //Fluid Force squared.
	double ueq[3] = { 0 };
	double weights[4] = { 0 };
	weights[0] = 0.296296296296296;// 8. / 27.;	
	weights[1] = 0.074074074074074;// 2. / 27.;
	weights[2] = 0.018518518518519;// 1. / 54.;
	weights[3] = 0.004629629629630;// 1. / 216.;
								   //double Fmaxsq = F_sum.maxsq();

								   /*
								   for (iz = 0; iz < Nz; iz++) {
								   for (iy = 0; iy < Ny; iy++) {
								   for (ix = 0; ix < Nx; ix++) {
								   erodelist(ix, iy, iz) = 0.;
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
								   erodelist(ix, iy, iz) = (-kappa_er*sqrt(FFsq - WDWsqsum)); //erode point. How to choose limit for erosion of a point? (m_star). Need to add code for F_sum. Only things that contribute should go into the sum.
								   }
								   else
								   erodelist(ix, iy, iz) = 0; //don't erode point
								   }
								   }
								   }
								   }
								   */


								   //erodes all points that's lost all their mass
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {

//				masschange(ix, iy, iz) += kappa_er*Delta_T*F_vdw(ix, iy, iz);
				if (solid_list(ix, iy, iz) == 0) { // if we're at a surface point
//					ero_reso_check(ix, iy, iz) += 1.;
					if (-masschange(ix, iy, iz) > masspernode) { // if more mass has eroded away than existing mass per point
//						if (ero_reso_check(ix, iy, iz) == 1.) { // if a point eroded away during one Delta_T window => too low resolution.
//							cout << "Resolution too low! Increase Delta_T, mass per node or Kappa.\n";
//							fprintf(errorfile, "%i", 1);
//						}
//						ero_reso_check(ix, iy, iz) = 0.;
						*ecp = 1;
						erodepoint(ix, iy, iz, masschange, solid_list, rho, u, nhat, ftemp, feq, e, F_vdw);
						/*masschange(ix, iy, iz) = 0.;
						solid_list(ix, iy, iz) = -1; //surface node becomes fluid node.
						rho(ix, iy, iz) = rho(ix + nhat(ix, iy, iz, 0), iy + nhat(ix, iy, iz, 1), iz + nhat(ix, iy, iz, 2)); //fluid node is initialized with same density as the interface node. This could prob be extrapolated.
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
						}*/
					}
				}

			}
		}
	}
	
	//erodes all points that become isolated. (Doesn't seem to work)
	/*
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				if (solid_list(ix, iy, iz) == 0) {	// surface node
					if (IsForeverAlone(ix, iy, iz, solid_list, e) == 1) // if surface node is isolated
						erodepoint(ix, iy, iz, masschange, solid_list, rho, nhat, f, e, F_vdw);
				}
			}
		}
	}
	*/


//	solid_list.printsolid_list(solfile);
	//masschange.clear();
	//F_sum.clear();
}
void erodepoint(int ix, int iy, int iz, density& masschange, Solid_list& solid_list, density& rho, velocity& u, Normalvector& nhat, direction_density& ftemp, EDF&feq, momentum_direction& e, vectorNcubed& F_vdw) {
	double weights[4] = { 0 };
	weights[0] = 0.296296296296296;// 8. / 27.;	
	weights[1] = 0.074074074074074;// 2. / 27.;
	weights[2] = 0.018518518518519;// 1. / 54.;
	weights[3] = 0.004629629629630;// 1. / 216.;
	int ixshift = 0;
	int iyshift = 0;
	int izshift = 0;
	int ashift = 0;
	int ixinterface = ix + nhat(ix,iy,iz,0);
	int iyinterface = iy + nhat(ix,iy,iz,1);
	int izinterface = iz + nhat(ix,iy,iz,2);
	double ueq[3] = { 0 };
	ueq[0] = u(ixinterface, iyinterface, izinterface, 0);
	ueq[1] = u(ixinterface, iyinterface, izinterface, 1);
	ueq[2] = u(ixinterface, iyinterface, izinterface, 2);
	masschange(ix, iy, iz) = 0.;
	solid_list(ix, iy, iz) = -1; //surface node becomes fluid node.
	rho(ix, iy, iz) = rho(ixinterface, iyinterface, izinterface); // fluid node is initialized with same density as the interface node. This could prob be extrapolated.
	//rho(ix, iy, iz) = 1.;																									 // check that rho != 0 at new points.
	for (int a = 0; a < 27; a++) { //Initialize new fluid point with same f as interface node. Also, check all nearby nodes. If it's a interior solid node, it becomes a surface.
		ftemp(ix, iy, iz, a) = ftemp(ixinterface, iyinterface, izinterface, a); // same f as interface node.
		feq(ix, iy,iz,a) = calcfeq(rho, e, ueq, ix, iy, iz, a, edfforcedir);
		ixshift = ix + e(a, 0);
		iyshift = iy + e(a, 1);
		izshift = iz + e(a, 2);
		//ashift = 2 * (13 - a);
		if (solid_list(ixshift, iyshift, izshift) == 1) { // if solid node, turn it to a surface node, stream, vdw force.
			solid_list(ixshift, iyshift, izshift) = 0;
			masschange(ixshift, iyshift, izshift) = 0.;
			//f(ixshift, iyshift, izshift, atemp + ashift) = f(ix + nhat(ix, iy, iz, 0), iy + nhat(ix, iy, iz, 1), iz + nhat(ix, iy, iz, 2), atemp);
			ftemp(ixshift, iyshift, izshift, a) = ftemp(ix,iy,iz,a);
			//rho(ixshift, iyshift, izshift) = 1.;
			computeVDWforce(ixshift, iyshift, izshift, e, solid_list, F_vdw);
		}
		else if (solid_list(ixshift, iyshift, izshift) == 0) {
			ftemp(ixshift, iyshift, izshift, a) = ftemp(ix,iy,iz, a);
			computeVDWforce(ixshift, iyshift, izshift, e, solid_list, F_vdw);
		}
	}
}
int IsForeverAlone(int ix, int iy, int iz, Solid_list& solid_list, momentum_direction& e) {
	//checks if solid (surface) node is isolated. 
	int ixshift;
	int iyshift;
	int izshift;
	int friends = 26;
	for (int a = 0; a < 27; a++) {		// look through all nearby nodes
		ixshift = ix + e(a, 0);
		iyshift = iy + e(a, 1);
		izshift = iz + e(a, 2);
		if (solid_list(ixshift, iyshift, izshift) == -1)
			friends--;
	}
	if (friends == 0)
		return 1;
	else
		return 0;
}
