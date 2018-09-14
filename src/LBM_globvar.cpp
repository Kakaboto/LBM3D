#pragma once
#ifndef LBM_GLOBVAR_H
#define LBM_GLOBVAR_H
#include "../LBM_globvar.h"
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

int Nx = 120;
int Ny = 50;
int Nz = 50;
int edfforcedir = 0; //0 = no gravity. 1 = x. 2 = y. 3 = z.
                                         // for the moment, only velocity BC in x direction.
string x0BC = "velocityBC";
string xNBC = "velocityBC";
string y0BC = "periodicBC";
string yNBC = "periodicBC";
string z0BC = "periodicBC";
string zNBC = "periodicBC";

double Bvel[3] = { 0.1, 0, 0 };

double latspace = 1.;
int tend = 2001;
//double erosiontol = 4e-8;
double kappa_er = 1.; // material property of solid. depends on toughness and density.
double WDWforce = 7.e-11; //Wan-Der-Waals force
int updatefreq = 700;
double c = 1.;
double tau = 0.502;
double mu = c*c*(tau - 0.5);
double umax_theo = 0.1;
double Re = 1;
double Lz = (Nz - 1)*latspace;
//double gg = mu*Re / (tau*L);
double gg = 0.001;
double F = mu * 10 / Lz;
int cellist[27] = { 3, 2, 3, 2, 1, 2, 3, 2, 3  ,  2, 1, 2, 1, 0, 1, 2, 1, 2  ,  3, 2, 3, 2, 1, 2, 3, 2, 3 };
int obchoice = 1;
#endif 
