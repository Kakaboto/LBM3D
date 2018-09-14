#pragma once
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "string.h"
#include "random"
#include "time.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <sstream>

using namespace std;
typedef vector<double> dvec;
//typedef Eigen::Matrix<double, 3, 3> Rotmat;

#define sq3inv 0.577350269189626
#define sq2inv 0.707106781186547
#define sq3 1.732050807568877
#define pi 3.1415926535897932384626433

extern int Nx;
extern int Ny;
extern int Nz;
extern int edfforcedir; //0 = no gravity. 1 = x. 2 = y. 3 = z.
                                         // for the moment, only velocity BC in x direction.
extern string x0BC;
extern string xNBC;
extern string y0BC;
extern string yNBC;
extern string z0BC;
extern string zNBC;

extern double Bvel[3];

extern double latspace;
extern int tend;
//double erosiontol = 4e-8;
extern double kappa_er; // material property of solid. depends on toughness and density.
extern double WDWforce; //Wan-Der-Waals force
extern int updatefreq;
extern double c;
extern double tau;
extern double mu;
extern double umax_theo;
extern double Re;
extern double Lz;
//double gg = mu*Re / (tau*L);
extern double gg;
extern double F;
extern int cellist[27];
extern int obchoice;
