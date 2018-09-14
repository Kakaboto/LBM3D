#ifndef FILE_globvar_SEEN
#define FILE_globvar_SEEN
/*---------------------------------------*/
#pragma once
// These header files can go to a common header */
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
//#include "LBM_utilityfunc.h"

using namespace std;
typedef vector<double> dvec;
//typedef Eigen::Matrix<double, 3, 3> Rotmat;

/* Following are definitions that all functions may need 
   they go to the same common hearder file */
#define sq3inv  0.577350269189626
#define sq2inv  0.707106781186547
#define sq3  1.732050807568877
#define pi  3.1415926535897932384626433

/* These are grid variables, they go to a structure called
   GrivVar*/
#define Nx  120
#define Ny  50
#define  Nz 50

string x0BC = "velocityBC";
string xNBC = "velocityBC";
string y0BC = "periodicBC";
string yNBC = "periodicBC";
string z0BC = "periodicBC";
string zNBC = "periodicBC";
double latspace = 1.;
/* These are parameters of the basic 
quantities needed by the LBM, they go to 
a structre LBMProp */

double Bvel[3] = { 0.1, 0, 0 };
double c = 1.;
double tau = 0.502;
double mu = c*c*(tau - 0.5);
double Re = 1;
double umax_theo = 0.1;
double Lz = (Nz - 1)*latspace;
double gg = 0.001;
double F = mu * 10 / Lz;
int edfforcedir = 0; //0 = no gravity. 1 = x. 2 = y. 3 = z.
					 // for the moment, only velocity BC in x direction.
int cellist[27] = { 3, 2, 3, 2, 1, 2, 3, 2, 3  ,  2, 1, 2, 1, 0, 1, 2, 1, 2  ,  3, 2, 3, 2, 1, 2, 3, 2, 3 };
//double gg = mu*Re / (tau*L);
/* These are timestepping parameter, they are part of
   TimeStep structure */
int tend = 2001;
//double erosiontol = 4e-8;
/* These are erosion properties; they 
   go to the structure EroProp */
double kappa_er = 1.; // material property of solid. depends on toughness and density.
double WDWforce = 7.e-11; //Wan-Der-Waals force
int updatefreq = 700;
#endif /* !FILE_globvar_SEEN */
