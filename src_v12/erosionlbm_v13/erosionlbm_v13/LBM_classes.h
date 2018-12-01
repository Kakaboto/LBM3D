#pragma once
#include "LBM_input.h"

class momentum_direction {
public:
	int **  e;
	momentum_direction() {
		e = new int*[27];
		//		for (auto &it : e) {
		//			it.resize(3);
		//		}
		for (int i = 0; i < 27; i++)
			e[i] = new int[3];
		int ix = 0;
		int iy = 0;
		int iz = 0;
		int a = 0;
		for (iz = -1; iz < 2; iz++) {
			for (iy = -1; iy < 2; iy++) {
				for (ix = -1; ix < 2; ix++) {
					e[a][0] = ix;
					e[a][1] = iy;
					e[a][2] = iz;
					//					cout << " " << e[a][0] << " " << e[a][1] << " " << e[a][2] << "\n";
					a++;
				}
			}
		}
	}

	int& operator()(int a, int i) {
		return e[a][i];
	}

private:
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;

};


class Grid {

public:
	dvec x;			// 0:Nx-1
	dvec y;			// 0:Nx-1
	dvec z;			// 0:Nx-1
	Grid() {
		x.resize(Nx);
		y.resize(Nx);
		z.resize(Nz);
		int ix = 0;
		int iy = 0;
		int iz = 0;
		//		latspace_x = 1. / (Nx - 1.);
		//		latspace_y = 1. / (Ny - 1.);
		//		latspace_z = 1. / (Nz - 1.);

		for (ix = 0; ix < Nx; ix++)
			x[ix] = xmin + latspace*ix;
		for (iy = 0; iy < Ny; iy++)
			y[iy] = ymin + latspace*iy;
		for (iz = 0; iz < Nz; iz++)
			z[iz] = zmin + latspace*iz;
		cout << "\n Grid created. \n";

	}
	void printgrid(FILE * gridfile) {
		int ix = 0;
		int iy = 0;
		int iz = 0;
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					fprintf(gridfile, "%e %e %e\n", x[ix], y[iy], z[iz]);
				}
			}
		}
		cout << "\n Printed grid to file for figures. \n";
	}
	~Grid();


private:
	double xmin = 0; double xmax = 1;
	double ymin = 0; double ymax = 1;
	double zmin = 0; double zmax = 1;
	double latspace_x = 0;
	double latspace_y = 0;
	double latspace_z = 0;
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;

};


class Solid_list {

public:
	vector<int> element;		// solid_list:	-1 = fluid point. 0 = surface solid point. 1 = interior solid point.
	Solid_list(int choice, Grid& grid, double rotation, momentum_direction& e, FILE * parfile);
	~Solid_list();
	int& operator()(int ix, int iy, int iz) {
		return element[(ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot];
	}

	double center[3];
	double center1[3];
	double center2[3];
	//double rotation = pi*0.5;

	void updatesurface(int ix, int iy, int iz, momentum_direction& e) {
		int a = 0;
		int n = 0;
		element[(ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot] = -1;
		for (a = 0; a < 27; a++) {
			n = (ix + 1 + e(a, 0)) + (iy + 1 + e(a, 1))*Nxtot + (iz + 1 + e(a, 2))*Nxtot*Nytot;
			if (element[n] == 1)
				element[n] = 0;
		}
	}

	void printsolid_list(FILE * solfile) {
		int ix = 0;
		int iy = 0;
		int iz = 0;
		int n = 0;
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					n = (ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot;
					fprintf(solfile, "%i ", element[n]);
				}
			}
		}
		fprintf(solfile, "\n");
		cout << "\n Printed grid to file. \n";
	}

	void clear_list() {
		int n = 0;
		for (int iz = -1; iz < Nz + 1; iz++) {
			for (int iy = -1; iy < Ny + 1; iy++) {
				for (int ix = -1; ix < Nx + 1; ix++) {
					n = (ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot;
					element[n] = 0;
				}
			}
		}
	}

private:
	//	double latspace = 1. / (Nx - 1); //assuming Nx = Ny = Nz
	double xmin = 0; double xmax = 1;
	double ymin = 0; double ymax = 1;
	double zmin = 0; double zmax = 1;
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;

};
// Unnessesary to fill list with 0:s at the start.



class Surface_fluid_list {
public:
	vector<int> element;
	Surface_fluid_list(Solid_list& solid_list, momentum_direction& e);
	~Surface_fluid_list();

	int& operator()(int ix, int iy, int iz) {
		return element[(ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot];
	}

private:
	//	double latspace = 1. / (Nx - 1); //assuming Nx = Ny = Nz
	double xmin = 0; double xmax = 1;
	double ymin = 0; double ymax = 1;
	double zmin = 0; double zmax = 1;
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;
};




class direction_density {
public:
	dvec element;		//size: 27*Ncubedtot
	direction_density() {
		element.resize(27 * Ncubedtot);
	}

	int index(int ix, int iy, int iz, int a) {
		int ind = (ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot + a*Ncubedtot;
		return ind;
	}

	double& operator()(int ix, int iy, int iz, int a) {
		return element[(ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot + a*Ncubedtot];
	}

	void clearvector() {
		// Fortsätt här!!
	}
private:
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;

};

class EDF {
public:
	dvec element;		//size: 27*Ncubed

	EDF() {
		element.resize(27 * Ncubed);
		fill(element.begin(), element.end(), 0);
	};

	int index(int ix, int iy, int iz, int a) {
		int ind = ix + iy*Nx + iz*Nx*Ny + a*Ncubed;
		return ind;
	}

	double& operator()(int ix, int iy, int iz, int a) {
		return element[ix + iy*Nx + iz*Nx*Ny + a*Ncubed];
	}
private:
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;
};

class velocity {
public:
	dvec element;		//size: 3*Ncubed


	velocity() {
		//element.resize(3 * Ncubed);
		element.resize(3 * Ncubed);
	};

	int index(int ix, int iy, int iz, int i) {
		int ind = ix + iy*Nx + iz*Nx*Ny + i;
		return ind;
	}

	double& operator()(int ix, int iy, int iz, int i) {
		return element[3 * (ix + iy*Nx + iz*Nx*Ny) + i];
	}

private:
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;
};

class density {
public:
	dvec element;		//size: Ncubed

	density() {
		element.resize(Ncubed);
		fill(element.begin(), element.end(), 1);
	};

	void clear() {
		fill(element.begin(), element.end(), 0.);
	}

	int index(int ix, int iy, int iz) {
		int ind = ix + iy*Nx + iz*Nx*Ny;
		return ind;
	}

	double& operator()(int ix, int iy, int iz) {
		return element[ix + iy*Nx + iz*Nx*Ny];
	}
	void operator()(int ix, int iy, int iz, double num) {
		element[ix + iy*Nx + iz*Nx*Ny] = num;
	}
private:
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;
};

class Stresstensor {
public:
	dvec element;
	Stresstensor() {
		element.resize(9 * Ncubed);
		fill(element.begin(), element.end(), 0);
	}
	~Stresstensor() {
		//		element.clear;
	}
	void clear() {
		fill(element.begin(), element.end(), 0);
	}

	double& operator()(int ix, int iy, int iz, int i, int j) {
		return element[(ix + iy*Nx + iz*Nx*Ny) * 9 + i * 3 + j];
	}
private:
	int Ncubed = Nx*Ny*Nz;

};

// same structure as velocity class
class Normalvector {
public:
	dvec element;
	Normalvector() {
		element.resize(3 * Ncubed);
		fill(element.begin(), element.end(), 0);
	}
	~Normalvector() {
		//		element.clear;
	}

	double& operator()(int ix, int iy, int iz, int i) {
		return element[(ix + iy*Nx + iz*Nx*Ny) * 3 + i];
	}
private:
	int Ncubed = Nx*Ny*Nz;
};
// same structure as velocity class
class Wall_force {
public:
	dvec element;		//size: Ncubed

	Wall_force() {
		element.resize(3 * Ncubed);
	};

	void clear() {
		fill(element.begin(), element.end(), 0.);
	}

	double& operator()(int ix, int iy, int iz, int i) {
		return element[(ix + iy*Nx + iz*Nx*Ny) * 3 + i];
	}

private:
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;
};

class vector3Ncubed {
public:
	dvec element;		//size: 3*Ncubed


	vector3Ncubed() {
		//element.resize(3 * Ncubed);
		element.resize(3 * Ncubed);
		fill(element.begin(), element.end(), 0.);
	};

	void clear() {
		fill(element.begin(), element.end(), 0.);
	}

	int index(int ix, int iy, int iz, int i) {
		int ind = 3 * (ix + iy*Nx + iz*Nx*Ny) + i;
		return ind;
	}

	double& operator()(int ix, int iy, int iz, int i) {
		return element[3 * (ix + iy*Nx + iz*Nx*Ny) + i];
	}

	double maxsq() {
		double temp[3] = { 0 };
		double max[3] = { 0 };
		int n = 0;
		for (int iz = 0; iz < Nz; iz++) {
			for (int iy = 0; iy < Ny; iy++) {
				for (int ix = 0; ix < Nx; ix++) {
					n = 3 * (ix + iy*Nx + iz*Nx*Ny);
					temp[0] = element[n];
					temp[1] = element[n + 1];
					temp[2] = element[n + 2];
					if (abs(temp[0]) > max[0])
						max[0] = temp[0];
					if (abs(temp[1]) > max[1])
						max[1] = temp[1];
					if (abs(temp[2]) > max[2])
						max[2] = temp[2];
				}
			}
		}
		return max[0] * max[0] + max[1] * max[1] + max[2] * max[2];
	}
private:
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;
};

class vectorNcubed {
public:
	dvec element;		//size: Ncubed

	vectorNcubed() {
		element.resize(Ncubed);
		fill(element.begin(), element.end(), 0.);
	};

	void clear() {
		fill(element.begin(), element.end(), 0.);
	}

	int index(int ix, int iy, int iz) {
		int ind = ix + iy*Nx + iz*Nx*Ny;
		return ind;
	}

	double& operator()(int ix, int iy, int iz) {
		return element[ix + iy*Nx + iz*Nx*Ny];
	}

private:
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;
};
