#pragma once

//#include <Eigen\Dense>
//#include "mpi.h"



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




// DATA LAYOUT FOR f AND ftemp:
// f: 
// (i,j,k):		(0,0,0):(Nx-1,Ny-1,Nz-1) - physical points. (-1,-1,-1), (Nx, Ny, Nz) - Boundary points ("ghost points").

// solid_list:	-1 = fluid point. 0 = surface solid point. 1 = interior solid point.

//Frågor:	Varför kan jag inte printa ut värden från solid_list å f när dom är const?

// Antar att Nx = Ny = Nz (isotropt grid).

// rho = 1 i början, även i solid points. Detta löser sig med att macrovariables sätter rätt saker till 0 i första iterationen. rho = 1
// för att velocity BC behöver initial density. 

// varför är c = 1/sqrt(3)?

// antar alltid att objekt ligger i mitten av lådan (sfären går att välja plats för)

// erosion: varje solid point antas bestå av 100st korn. Each erosion step determines how many grains get eroded away (if  all 100 goes away then we say the point has eroded).



using namespace std;
typedef vector<double> dvec;


//typedef Eigen::Matrix<double, 3, 3> Rotmat;

extern double sq3inv;
extern double sq2inv;
extern double sq3;
extern double sq2;
extern double pi;

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

//double erosiontol = 4e-8;
extern double masspernode; // Amount of mass per solid node. units = kg?
extern double kappa_er; // material property of solid. depends on toughness and density.
extern double VDW_0; //Wan-Der-Waals force
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
extern double Delta_T; //Need Delta_T to be bigger than characteristic time. Hence, the factor 10 in front of L/U.
extern double dt;
extern int obchoice; 
