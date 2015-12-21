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
	 * 0 = simple geometry with variable Re number
	 */
	int geom_type(-1);
	double alpha(M_PI/4.0);
	double re(10000);
	int id(0);
	
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
				if (geom_type == 4) {
					// read alpha (used in the Karman vortex street)
					alpha = atof(argv[i]);
				} else if (geom_type == 0) {
					re = atof(argv[i]);
				}
				break;
			case 4:
				id = atoi(argv[i]);
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
		case 0:
			geom_gen.simpleGeom(re);
			break;
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
		case 6:
			geom_gen.testCase2();
			break;
		case 7:
			geom_gen.test_twoCellCriterion();
			break;
		case 8:
			geom_gen.test_twoCellCriterion2();
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
	if (geom_type == 0) {
		std::string file1("../uq_data/uq_parameter_" + std::to_string(id) + ".params");
		geom_gen.writeToFile_simple("../uq_data/uq_geometry.geom", file1.c_str());
	} else {
		geom_gen.writeToFile("../data/complex_default.geom", "../data/assignment03.params");
	}

}
