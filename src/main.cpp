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
#include "typedef.hpp"
#include "communicator.hpp"
#include "compute.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "visu.hpp"
#include "vtk.hpp"

#include <iostream>
#include <sys/stat.h>

#ifdef __linux__
	#include <sys/time.h>
//#elif _WIN32 // windows 32- and 64-bit
//	#include <Windows.h>
#endif

#define VERBOSE false

int main(int argc, char **argv) {

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

  // Create communicator
  Communicator comm(&argc, &argv); // TODO: handle arguments

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
  Parameter param;
  Geometry geom(&comm);

	// do domain decomposition on master and send information to all other processes
	multi_index_t tdim;
	int** rankDistri;
	multi_index_t** localSizes;
	if (comm.getRank() == 0) {
		std::cout << "Demain decomposition...\n" << std::flush;
		// horizontal decomposition
		index_t np = comm.getSize();
		tdim[0] = np;
		tdim[1] = 1;
		// allocate
		rankDistri = new int*[tdim[0]];
		for (int i=0; i<tdim[0]; i++){
			rankDistri[i] = new int[tdim[1]];
		}
		localSizes = new multi_index_t*[tdim[0]];
		for (int i=0; i<tdim[0]; i++){
			localSizes[i] = new multi_index_t[tdim[1]];
		}
		// write values
		for (int i=0; i<np; i++){
			rankDistri[i][0] = i;
			localSizes[i][0] = multi_index_t(geom.TotalSize()[0]/real_t(np), geom.TotalSize()[1]);
		}
		// send information
		//std::cout << "Broadcasting domain decomposition information...\n" << std::flush;
		int sendBuf(tdim[0]);
		MPI_Bcast(&sendBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
		sendBuf = tdim[1];
		MPI_Bcast(&sendBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);

		//std::cout << "First broadcasting completed.\n" << std::flush;

		for (int i=0; i<tdim[0]; i++){
			MPI_Bcast(rankDistri[i], tdim[1], MPI_INT, 0, MPI_COMM_WORLD);
		}

		//std::cout << "Second broadcasting completed.\n" << std::flush;

		for (int i=0; i<tdim[0]; i++){
			for (int j=0; j<tdim[1]; j++) {
				sendBuf = localSizes[i][j][0];
				MPI_Bcast(&sendBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
				sendBuf = localSizes[i][j][1];
				MPI_Bcast(&sendBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		//std::cout << "All broadcasting completed.\n" << std::flush;
	} else {
		// receive information
		int recBuf(0);
		MPI_Bcast(&recBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
		tdim[0] = recBuf;
		MPI_Bcast(&recBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
		tdim[1] = recBuf;

		//std::cout << "Recieved tdim: " << tdim[0] << ", " << tdim[1] << "\n" << std::flush;

		// allocate
		rankDistri = new int*[tdim[0]];
		for (int i=0; i<tdim[0]; i++){
			rankDistri[i] = new int[tdim[1]];
		}
		localSizes = new multi_index_t*[tdim[0]];
		for (int i=0; i<tdim[0]; i++){
			localSizes[i] = new multi_index_t[tdim[1]];
		}

		for (int i=0; i<tdim[0]; i++){
			MPI_Bcast(rankDistri[i], tdim[1], MPI_INT, 0, MPI_COMM_WORLD);
		}
		for (int i=0; i<tdim[0]; i++){
			for (int j=0; j<tdim[1]; j++) {
				MPI_Bcast(&recBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
				localSizes[i][j][0] = recBuf;
				MPI_Bcast(&recBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
				localSizes[i][j][1] = recBuf;
				//std::cout << "Recieved localSizes[i][j]: " << localSizes[i][j][0] << ", " << localSizes[i][j][1] << "\n" << std::flush;
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	// write to communicator
	comm.setProcDistribution(rankDistri, tdim, localSizes);

	// update geometry values
	geom.update_values();

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
  Compute comp(&geom, &param, &comm);
	if (VERBOSE) std::cout << "Done.\n" << std::flush;

  if (comm.getRank() == 0) { // TODO: do this on every process? what if the processes run on different nodes?
    // check if folder "VTK" exists
    struct stat info;

    if (stat("VTK", &info) != 0) {
      system("mkdir VTK");
    }
  }

// Create and initialize the visualization
#ifdef USE_DEBUG_VISU
	if (VERBOSE) std::cout << "Initializing the visualization..." << std::flush;
  Renderer visu(geom.Length(), geom.Mesh());
  visu.Init(800 / comm.ThreadDim()[0], 800 / comm.ThreadDim()[1],
            comm.getRank() + 1);
	if (VERBOSE) std::cout << "Done.\n" << std::flush;
#endif // USE_DEBUG_VISU

  // Create a VTK generator;
	if (VERBOSE) std::cout << "Creating the VTK generator..." << std::flush;
  // use offset as the domain shift
  multi_real_t offset;
  offset[0] = comm.ThreadIdx()[0] * (geom.Mesh()[0] * (double)(geom.Size()[0] - 2));
  offset[1] = comm.ThreadIdx()[1] * (geom.Mesh()[1] * (double)(geom.Size()[1] - 2));
  VTK vtk(geom.Mesh(), geom.Size(), geom.TotalSize(), offset, comm.getRank(),
          comm.getSize(), comm.ThreadDim());
	if (VERBOSE) std::cout << "Done.\n" << std::flush;

#ifdef USE_DEBUG_VISU
  const Grid *visugrid;

  visugrid = comp.GetVelocity();
#endif // USE_DEBUG_VISU

  // Run the time steps until the end is reached
  while (comp.GetTime() < param.Tend()) {
#ifdef USE_DEBUG_VISU
    // Render and check if window is closed
    switch (visu.Render(visugrid)) {
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
    default:
      break;
    };
	visugrid->CheckNaN(); // check for NaNs
#endif // USE_DEBUG_VISU

    // Create VTK Files in the folder VTK
	if (VERBOSE) std::cout << "Creating VTK files..." << std::flush;
    // Note that when using VTK module as it is you first have to write cell
    // information, then call SwitchToPointData(), and then write point data.
    vtk.Init("VTK/field");
    vtk.AddRank();
    vtk.AddCellField("Cell Velocity", comp.GetU(), comp.GetV());
    vtk.SwitchToPointData();
    vtk.AddPointField("Velocity", comp.GetU(), comp.GetV());
    vtk.AddPointScalar("Pressure", comp.GetP());
    vtk.Finish();
	if (VERBOSE) std::cout << "Done.\n" << std::flush;

    // Run a few steps
	if (VERBOSE) std::cout << "Running a few timesteps...\n" << std::flush;
    for (uint32_t i = 0; i < 9; ++i)
      comp.TimeStep(false,VERBOSE);
    bool printOnlyOnMaster = !comm.getRank();
    comp.TimeStep(printOnlyOnMaster,VERBOSE);
  }

	std::cout << "The program was executed sucessfully!\n";

	// runtime measurement
	if (comm.getRank() == 0) { // do this only on the master
#ifdef __linux__
		gettimeofday(&tv, NULL);
		milliseconds_end = tv.tv_sec * 1000 + tv.tv_usec / 1000;
		diff_sec = (milliseconds_end - milliseconds_begin) / 1000.0;
		std::cout << "Total elapsed time: " << diff_sec << " seconds.\n";
//#elif _WIN32 // windows 32- and 64-bit
//	
#endif
	}

	std::cout << "Exiting...\n" << std::flush;

  return 0;
}
