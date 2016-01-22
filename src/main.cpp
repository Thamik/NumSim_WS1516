/*
 * Copyright (C) 2015   Malte Brunn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//------------------------------------------------------------------------------

//#define VTK_OUTPUT
#define OUTPUT_MAIN

#include "typedef.hpp"
#include "communicator.hpp"
#include "geometry.hpp"
#include "solver.hpp"
#include "visu.hpp"

#ifdef VTK_OUTPUT
#include "vtk.hpp"
#endif

#include <iostream>
#include <sys/stat.h>
#include <stdlib.h>
#include <fstream>

#include <math.h>
#define _USE_MATH_DEFINES

#ifdef __linux__
	#include <sys/time.h>
//#elif _WIN32 // windows 32- and 64-bit
//	#include <Windows.h>
#endif

int main(int argc, char **argv) {

	// read the command line arguments
	std::string param_file("");
	std::string geom_file("");
	std::string uq_filename("");
#ifdef UQ_RUN
	int sim_id(0);
#endif
	for (int i=0; i<argc; i++){
		switch (i){
			case 0:
				// this is the path to the executable
				break;
			case 1:
				// this should be the filename of the geometry file
				geom_file.assign(argv[i]);
				break;
			default:
				// too much parameters, don't know what to do here
				std::cout << "Warning: too much command line arguments!\n";
				break;
		}
	}

	// Create communicator
	Communicator comm(&argc, &argv); // the arguments are being handed over to the MPI_Init call	

	Geometry geom(&comm);

	// load data from files, if filenames given in the command line arguments
	if (geom_file.compare("") != 0){
		geom.load_domain_partitioning(geom_file.c_str());
	} else {
		geom.load_domain_partitioning(nullptr);
	}

	Grid* _p = new Grid(&geom, multi_real_t(-0.5,-0.5));
	Grid* _rhs = new Grid(&geom, multi_real_t(-0.5,-0.5));
	MGSolver* _solver = new MGSolver(pow(10,-4),1000);
	real_t residual(0.0);
	/*real_t omega = 1.7;
	RedOrBlackSOR* _solver = new RedOrBlackSOR(&geom, omega);
	real_t residual;
	index_t totalIteration(0);*/

	// measure runtime
	timeval tv;
	long int milliseconds_begin(0), milliseconds_end(0);
	double diff_sec(0.0);
	if (comm.getRank() == 0) { // do this only on the master process
#ifdef __linux__
		gettimeofday(&tv, NULL);
		milliseconds_begin = tv.tv_sec * 1000 + tv.tv_usec / 1000;
//#elif _WIN32 // windows 32- and 64-bit
//
#else
		std::cout << "System unknown, no runtime measurement!\n" << std::flush;
#endif
	}
	
	for (int ii=0; ii<10; ii++) {
		_p->Initialize(1.0);
		_rhs->Initialize(0.0);
		//MGSolver* _solver = new MGSolver(pow(10,-4), 1000);
		const MGInfoHandle info = _solver->Solve(_p, _rhs);
		residual = info.getResidual();
		/*residual = (pow(10,-4) + 1.0);
		index_t iteration(0);
		while (true){
			// one solver cycle is done here, alternating red or black
			if ((iteration % 2) == (comm.EvenOdd() ? 1 : 0)){
				residual = _solver->RedCycle(_p, _rhs);
			} else {
				residual = _solver->BlackCycle(_p, _rhs);
			}

			iteration++;

			// check for edge condition
			if (iteration > 1000000*2){
				break;
			} else if (residual < pow(10,-4)){
				break;
			}
		}
		totalIteration += iteration;*/
	}
	// _solver_converging = info.getConverged();

	// runtime measurement
	if (comm.getRank() == 0) { // do this only on the master

#ifdef OUTPUT_MAIN
		std::cout << "The program was executed sucessfully!\n";

#ifdef __linux__
		gettimeofday(&tv, NULL);
		milliseconds_end = tv.tv_sec * 1000 + tv.tv_usec / 1000;
		diff_sec = (milliseconds_end - milliseconds_begin) / 1000.0;
		std::cout << "Total elapsed time (for 10 Simulations): " << diff_sec << " seconds.\n";
//#elif _WIN32 // windows 32- and 64-bit
//	
#endif
		std::cout << "Exiting...\n" << std::flush;
#endif	

		// delete variables
		delete _p;
		delete _rhs;
		delete _solver;

		//std::cout << "Iteration number: " << totalIteration << std::endl;
		std::cout << "Solver converged with residuum: " << residual << std::endl;
	}


	return 0;
}
