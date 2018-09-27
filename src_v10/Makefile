default:
	g++ -std=c++11 -c LBM_input.cpp
	g++ -std=c++11 -c LBM_classes.cpp
	g++ -std=c++11 -c LBM_utility_func.cpp
	g++ -std=c++11 -c LBM_fluidsteps_func.cpp
	g++ -std=c++11 -c LBM_force_erosion_func.cpp
	g++ -std=c++11 -c erosionlbm_v11.cpp
	g++ LBM_input.o LBM_classes.o LBM_utility_func.o LBM_fluidsteps_func.o LBM_force_erosion_func.o erosionlbm_v11.o -o ../Re01/ero_Re01.ext

clean:
	rm -f *.out
