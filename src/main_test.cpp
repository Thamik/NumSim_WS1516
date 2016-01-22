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

#define VTK_OUTPUT
#define OUTPUT_MAIN

#include "typedef.hpp"
#include "communicator.hpp"
#include "compute.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "visu.hpp"

#ifdef VTK_OUTPUT
#include "vtk.hpp"
#endif

#include <iostream>
#include <sys/stat.h>
#include <stdlib.h>
#include <fstream>

#ifdef __linux__
	#include <sys/time.h>
//#elif _WIN32 // windows 32- and 64-bit
//	#include <Windows.h>
#endif

int main(int argc, char **argv) {

	// Create communicator
	Communicator comm(&argc, &argv); // the arguments are being handed over to the MPI_Init call

	// print welcome screen
	if (!comm.getRank()){
		std::ifstream input;
		input.open("numsim_logo.txt");
		if (input.is_open()){
			std::string line;
			while(std::getline(input,line)) {
				std::cout << line << '\n';
			}
		}
		input.close();
	}

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
				// this should be the filename of the parameter file
				param_file.assign(argv[i]);
				break;
			case 2:
				// this should be the filename of the geometry file
				geom_file.assign(argv[i]);
				break;
#ifdef UQ_RUN
			case 3:
				// this should be the simulation id
				sim_id = atoi(argv[i]);
				break;
#endif
			default:
				// too much parameters, don't know what to do here
				std::cout << "Warning: too much command line arguments!\n";
				break;
		}
	}

#ifdef UQ_RUN
	uq_filename += "VTK/uq_data_" + std::to_string(sim_id) + ".uqdat";
#endif

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

	// Create parameter and geometry instances with default values
	Parameter param(&comm);
	Geometry geom(&comm);

	// load data from files, if filenames given in the command line arguments
	if (param_file.compare("") != 0){
		// the string is not empty, try to load the data
		param.Load(param_file.c_str());
	}
	if (geom_file.compare("") != 0){
		geom.load_domain_partitioning(geom_file.c_str());
	} else {
		geom.load_domain_partitioning(nullptr);
	}

	// Create the fluid solver
#ifdef UQ_RUN
	Compute comp(&geom, &param, &comm, uq_filename.c_str(), sim_id);
#else
	Compute comp(&geom, &param, &comm);
#endif

	if (comm.getRank() == 0) { // TODO: do this on every process? what if the processes run on different nodes?
		// check if folder "VTK" exists
		struct stat info;

		if (stat("VTK", &info) != 0) {
		int res_sys = system("mkdir VTK");
		if (!res_sys) std::cout << "Warning: mkdir VTK may have not worked!\n" << std::flush;
		}
	}

	// Create and initialize the visualization
#ifdef USE_DEBUG_VISU
	Renderer visu(geom.Length(), geom.Mesh(), &geom);

	index_t winMaxX = 1200;
	index_t winMaxY = 650;
	index_t winSizeX, winSizeY;
	if (geom.TotalLength()[0]/double(winMaxX) >= geom.TotalLength()[1]/double(winMaxY)){
		winSizeX = winMaxX;
		winSizeY = int(round(winSizeX * geom.TotalLength()[1] / geom.TotalLength()[0]));
	} else {
		winSizeY = winMaxY;
		winSizeX = int(round(winSizeY * geom.TotalLength()[0] / geom.TotalLength()[1]));
	}
	index_t sizeX = int(round(winSizeX * geom.Size()[0] / double(geom.TotalSize()[0])));
	index_t sizeY = int(round(winSizeY * geom.Size()[1] / double(geom.TotalSize()[1])));
	visu.Init(sizeX, sizeY, comm.getRank() + 1, comm.ThreadIdx(), comm.ThreadDim());
#endif // USE_DEBUG_VISU

#ifdef VTK_OUTPUT
	// Create a VTK generator;
	// use offset as the domain shift
	multi_real_t offset;

	offset[0] = geom.Mesh()[0] * double(geom.TotalOffset()[0]);
	offset[1] = geom.Mesh()[1] * double(geom.TotalOffset()[1]);

	VTK vtk(geom.Mesh(), geom.Size(), geom.TotalSize(), offset, comm.getRank(), comm.getSize(), comm.ThreadDim());
#endif

#ifdef USE_DEBUG_VISU
	const Grid *visugrid;

	visugrid = comp.GetVelocity();
#endif // USE_DEBUG_VISU

	// prepare console output
#ifdef OUTPUT_MAIN
	if(!comm.getRank()){
		std::cout << "\n\n\n\n\n" << "\n\n\n" << std::flush;
#ifdef USE_PARTICLES
		std::cout << "\n";
#endif
#ifndef GO_FAST
		std::cout << "\n\n";
#endif
	}
#endif

	// Run the time steps until the end is reached
	while (comp.GetTime() < param.Tend()) {
#ifdef USE_DEBUG_VISU
		// Render and check if window is closed

		real_t min_visugrid = visugrid->TotalInnerMin();
		real_t max_visugrid = visugrid->TotalInnerMax();
		switch (visu.Render(visugrid, min_visugrid, max_visugrid)) {
			case -1:
				return -1;
			case 0:
				visugrid = comp.GetVelocity();
				break;
			case 1:
				visugrid = comp.GetU();
				break;
			case 2:
				visugrid = comp.GetV();
				break;
			case 3:
				visugrid = comp.GetP();
				break;
			case 4:
				visugrid = comp.GetT();
				break;
			case 5:
				visugrid = comp.GetVorticity();
				break;
			default:
				break;
		};
#endif // USE_DEBUG_VISU

#ifdef VTK_OUTPUT
		// Create VTK Files in the folder VTK

		// Note that when using VTK module as it is you first have to write cell
		// information, then call SwitchToPointData(), and then write point data.

		vtk.Init("VTK/field", &geom);
		vtk.AddRank();
		vtk.AddCellField("Cell Velocity", comp.GetU(), comp.GetV());
		vtk.SwitchToPointData();
		vtk.AddPointField("Velocity", comp.GetU(), comp.GetV());
		vtk.AddPointScalar("Pressure", comp.GetP());
		vtk.AddPointScalar("Vorticity", comp.GetVorticity());
		vtk.AddPointScalar("Stream Function", comp.GetStream());
		vtk.AddPointScalar("Temperature", comp.GetT());
		vtk.Finish();
#endif

		// Run a few steps
		/*for (int i=0; i<25; i++){
			comp.TimeStep(false);
#ifdef UQ_RUN
			if (comm.getRank() == 0) comp.writeUQFile();
#endif
		}*/
		comp.TimeStep(true); // print info

	}

	// runtime measurement
	if (comm.getRank() == 0) { // do this only on the master

#ifdef OUTPUT_MAIN
		std::cout << "The program was executed sucessfully!\n";

#ifdef __linux__
		gettimeofday(&tv, NULL);
		milliseconds_end = tv.tv_sec * 1000 + tv.tv_usec / 1000;
		diff_sec = (milliseconds_end - milliseconds_begin) / 1000.0;
		std::cout << "Total elapsed time: " << diff_sec << " seconds.\n";
//#elif _WIN32 // windows 32- and 64-bit
//	
#endif

		std::cout << "Exiting...\n" << std::flush;
#endif
	}

	return 0;
}
