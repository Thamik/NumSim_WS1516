/*
 *  Copyright (C) 2015   Malte Brunn
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "typedef.hpp"
#include "compute.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "visu.hpp"
#include "vtk.hpp"

#include <iostream>

#ifdef __linux__
	#include <sys/time.h>
//#elif _WIN32 // windows 32- and 64-bit
//	#include <Windows.h>
#endif

#define VERBOSE false

int main(int argc, char **argv) {

	// measure runtime
#ifdef __linux__
	timeval tv;
	gettimeofday(&tv, NULL);
	long int milliseconds_begin = tv.tv_sec * 1000 + tv.tv_usec / 1000;
//#elif _WIN32 // windows 32- and 64-bit
//
#else
	std::cout << "System unknown, no runtime measuring!\n" << std::flush;
#endif

	// read the command line arguments
	std::string param_file("");
	std::string geom_file("");
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
			default:
				// too much parameters, don't know what to do here
				std::cout << "Warning: too much command line arguments!\n";
				break;
		}
	}

  // Create parameter and geometry instances with default values
  Parameter param;
  Geometry geom;

	// load data from files, if filenames given in the command line arguments
	if (param_file.compare("") != 0){
		// the string is not empty, try to load the data
		param.Load(param_file.c_str(), VERBOSE);
	}
	if (geom_file.compare("") != 0){
		// the string is not empty, try to load the data
		geom.Load(geom_file.c_str(), VERBOSE);
	}

  // Create the fluid solver
	if (VERBOSE) std::cout << "Creating the fluid solver..." << std::flush;
  Compute comp(&geom, &param);
	if (VERBOSE) std::cout << "Done.\n" << std::flush;

  // Create and initialize the visualization
	if (VERBOSE) std::cout << "Initializing the visualization..." << std::flush;
  Renderer visu(geom.Length(), geom.Mesh());
  visu.Init(800, 800);
	if (VERBOSE) std::cout << "Done.\n" << std::flush;

  // Create a VTK generator
	if (VERBOSE) std::cout << "Creating the VTK generator..." << std::flush;
  VTK vtk(geom.Mesh(), geom.Size());
	if (VERBOSE) std::cout << "Done.\n" << std::flush;

  const Grid *visugrid;
  bool run = true;

  visugrid = comp.GetP();

  // Run the time steps until the end is reached
  while (comp.GetTime() < param.Tend() && run) {
    // Render and check if window is closed
    switch (visu.Render(visugrid)) {
    case -1:
      run = false;
      break;
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
    default:
      break;
    };
	visugrid->CheckNaN();

    // Create a VTK File in the folder VTK (must exist)
	if (VERBOSE) std::cout << "Creating a VTK file...";
    vtk.Init("VTK/field");
    vtk.AddField("Velocity", comp.GetU(), comp.GetV());
    vtk.AddScalar("Pressure", comp.GetP());
    vtk.Finish();
	if (VERBOSE) std::cout << "Done.\n";

    // Run a few steps
	if (VERBOSE) std::cout << "Running a few timesteps...\n" << std::flush;
    for (uint32_t i = 0; i < 9; ++i){
      comp.TimeStep(false,VERBOSE);
    }
    comp.TimeStep(true,VERBOSE);
  }

	std::cout << "The program was executed sucessfully!\n";

#ifdef __linux__
	gettimeofday(&tv, NULL);
	long int milliseconds_end = tv.tv_sec * 1000 + tv.tv_usec / 1000;
	double diff_sec = (milliseconds_end - milliseconds_begin) / 1000.0;
	std::cout << "Total elapsed time: " << diff_sec << " seconds.\n";
//#elif _WIN32 // windows 32- and 64-bit
//	
#endif

	std::cout << "Exiting...\n" << std::flush;

  return 0;
}
