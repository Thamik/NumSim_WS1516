#include <iostream>
#include <stdlib.h>     /* atof */

#define _USE_MATH_DEFINES // to get pi via M_PI
#include <math.h>

#include "geom_gen.hpp"

int main(int argc, char** argv){

	// read command line arguments
	int totalSize(-1);
	/* geom_type:
	 * 1 = driven cavity
	 * 2 = pipe flow
	 * 3 = flow over a step
	 * 4 = karman vortex street
	 */
	int geom_type(-1);
	double alpha(M_PI/4.0);
	
	for (int i=0; i<argc; i++){
		switch (i){
			case 0:
				// this argument is the name of the executable
				break;
			case 1:
				totalSize = atoi(argv[i]);
				break;
			case 2:
				geom_type = atoi(argv[i]);
				break;
			case 3:
				// read alpha (used in the Karman vortex street)
				alpha = atof(argv[i]);
				break;
			default:
				std::cout << "Warning: too much command line arguments!\n" << std::flush;
		}
	}

	GeometryGenerator geom_gen;

	if (totalSize > 0){
		// set the read size
		geom_gen.setTotalSize(totalSize);
	} else {
		// set the default size
		geom_gen.setTotalSize(10000);
	}

	switch (geom_type){
		case 1:
			geom_gen.drivenCavity();
			break;
		case 2:
			geom_gen.pipeFlow();
			break;
		case 3:
			geom_gen.flowOverAStep();
			break;
		case 4:
			geom_gen.karmanVortexStreet(alpha);
			break;
		case 5:
			geom_gen.testCase1();
			break;
		case -1:
			// default case: no console input
			geom_gen.drivenCavity();
			break;
		default:
			std::cout << "Warning: Unknown geometry type! Using default..." << std::endl;
			geom_gen.drivenCavity();
	}

	geom_gen.print();
	geom_gen.writeToFile("../data/complex_default.geom");

}
