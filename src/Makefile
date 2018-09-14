default:
	g++ -std=c++11 -c LBM_globvar.cpp
	g++ -std=c++11 -c LBM_classes.cpp 
	g++ -std=c++11 -c LBM_utilityfunc.cpp 
	g++ -std=c++11 -c LBM_stepsfunc.cpp 
	g++ -std=c++11 -c LBM_force_erofunc.cpp 
	g++ -std=c++11 -c erosion.cpp 
	g++ LBM_globvar.o LBM_classes.o LBM_utilityfunc.o LBM_stepsfunc.o LBM_force_erofunc.o erosion.o -o ero.ext
clean:
	rm -f *.out
