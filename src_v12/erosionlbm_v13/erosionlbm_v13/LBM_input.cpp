#include "LBM_input.h"


//#include <Eigen\Dense>
//#include "mpi.h"







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



//typedef Eigen::Matrix<double, 3, 3> Rotmat;

double sq3inv = 0.577350269189626;
double sq2inv = 0.707106781186547;
double sq3 = 1.732050807568877;
double sq2 = 1.414213562373095;
double pi = 3.1415926535897932384626433;

int Nx = 120;
int Ny = 60;
int Nz = 60;
int edfforcedir = 0; //0 = no gravity. 1 = x. 2 = y. 3 = z.
					 // for the moment, only velocity BC in x direction.
string x0BC = "velocityBC";
string xNBC = "velocityBC";
string y0BC = "periodicBC";
string yNBC = "periodicBC";
string z0BC = "periodicBC";
string zNBC = "periodicBC";

double Bvel[3] = { 0.02, 0, 0 }; //inlet velocity

double latspace = 1.;

// for v = 0.01, tau = 0.5015, mpn = 5e-4 seems ok
// for v = 0.1: tau = 0.515, mpn = 0.1, VDW0 = 2.5*1e-5, kappa_er = 1 seems ok.
double masspernode = 0.1; // Amount of mass per solid node. units = kg?
double kappa_er = 1.; // material property of solid. depends on toughness and density. Lower value, higher toughness. Higher value, lower toughness.
double VDW_0 = 3*1e-5; //Wan-Der-Waals force
//int updatefreq = 700;
double c = 1.;
double tau = 1.;
double mu = 1./3.*(tau - 0.5);
double umax_theo = 0.1;
double Lz = (Nz - 1)*latspace;
//double gg = mu*Re / (tau*L);
double gg = 0.001;
double F = mu * 10 / Lz;
int cellist[27] = { 3, 2, 3, 2, 1, 2, 3, 2, 3  ,  2, 1, 2, 1, 0, 1, 2, 1, 2  ,  3, 2, 3, 2, 1, 2, 3, 2, 3 };
double dt = 1.;
int obchoice = 1;
double sphere_radius = 15.;
double Delta_T = 3.*2.*sphere_radius/Bvel[1]; //Need Delta_T to be bigger than characteristic time. Hence, the factor in front of L/U.
double Re = Bvel[0] * 2. * sphere_radius / mu;
//double Re = Bvel[0] * (double)Nz*2. / (3.*mu);
int tend = 1600;
