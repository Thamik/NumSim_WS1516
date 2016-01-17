#include <iostream>
#include <stdlib.h>     /* atof */

#define OUTPUT

#define _USE_MATH_DEFINES // to get pi via M_PI
#include <math.h>

#include <string>
#include <algorithm>

#include "geom_gen.hpp"

void frontend_interactive()
{
	GeometryGenerator geom_gen;
	int id(0);

	std::cout << "Welcome to NumSim Geometry Generator!\n";

	int geom_type(-1);
	std::cout << "Please choose the geometry type:\n";
	std::cout << "0 : simple geometry with variable Reynold number (uq run)\n";
	std::cout << "1 : driven cavity\n";
	std::cout << "2 : pipe flow\n";
	std::cout << "3 : flow over a step\n";
	std::cout << "4 : karman vortex street\n";
	std::cout << "5 : test case 1\n";
	std::cout << "6 : test case 2\n";
	std::cout << "7 : two cell criterion\n";
	std::cout << "8 : two cell criterion 2\n";
	while (true){
		std::cout << ">> ";
		std::cin >> geom_type;
		if (geom_type < 0 || geom_type > 8){
			std::cout << "Unknown geometry type, please type another number!\n";
		} else {
			break;
		}
	}

	int totalSize(-1);
	std::cout << "Please type the total number of cells:\n";
	while (true){
		std::cout << ">> ";
		std::cin >> totalSize;
		if (totalSize <= 0){
			std::cout << "The total number of cells can't be negative, please type another number!\n";
		} else {
			break;
		}
	}

	std::string inp;
	std::cout << "Would you like to force the number of cells in each dimension to be a power of 2? (y/n)\n";
	while (true){
		std::cout << ">> ";
		std::cin >> inp;
		if (inp.compare("y") == 0){
			geom_gen.forcePowerOfTwo();
			break;
		} else if (inp.compare("n") == 0){
			geom_gen.doNotForcePowerOfTwo();
			break;
		} else {
			std::cout << "Please type \"y\" or \"n\"!\n";
		}
	}

	geom_gen.setTotalSize(totalSize);

	if (geom_type == 0){
		std::cout << "Please type the Reynolds number:\n";
		std::cout << ">> ";
		double re(0.0);
		std::cin >> re;
		geom_gen.drivenCavity();
		geom_gen.setRe(re);

		std::cout << "Please type the sample id:\n";
		std::cout << ">> ";
		std::cin >> id;
	} else if (geom_type == 1){
		geom_gen.drivenCavity();
	} else if (geom_type == 2){
		geom_gen.pipeFlow();
	} else if (geom_type == 3){
		geom_gen.flowOverAStep();
	} else if (geom_type == 4){
		std::cout << "Please type the angle of the obstacle (in degree):\n";
		std::cout << ">> ";
		double alpha(0.0);
		std::cin >> alpha;
		alpha = alpha/360.0 * 2*M_PI; // transform to radiants
		geom_gen.karmanVortexStreet(alpha);
	} else if (geom_type == 5){
		geom_gen.testCase1();
	} else if (geom_type == 6){
		geom_gen.testCase2();
	} else if (geom_type == 7){
		geom_gen.test_twoCellCriterion();
	} else if (geom_type == 8){
		geom_gen.test_twoCellCriterion2();
	} else {
		// this should not happen, invalid geom_type
		throw std::runtime_error("Fatal error: geometry type invalid!");
	}

#ifdef OUTPUT
	geom_gen.print();
#endif

	if (geom_type == 0){
		std::string filenameParam("../uq_data/uq_parameter_" + std::to_string(id) + ".param");
		geom_gen.setFilenames("../uq_data/uq_geometry.geom", filenameParam.c_str());
	}

	geom_gen.writeToFile();
}

void frontend_commandlineargs(int argc, char** argv){
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
	double re(1500.0);
	int id(0);
	bool powerTwo(false);
	
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
				} else {
					powerTwo = true;
					std::cout << "Forcing pows of 2 enabled!" << std::endl;
				}
				break;
			case 4:
				id = atoi(argv[i]);
				break;
			default:
				std::cout << "Warning: too much command line arguments!\n" << std::flush;
		}
	}

	GeometryGenerator geom_gen;

	if (powerTwo){
		geom_gen.forcePowerOfTwo();
	} else {
		geom_gen.doNotForcePowerOfTwo();
	}

	if (totalSize > 0){
		// set the read size
		geom_gen.setTotalSize(totalSize);
	} else {
		// set the default size
		geom_gen.setTotalSize(10000);
	}

	switch (geom_type){
		case 0:
			geom_gen.drivenCavity();
			geom_gen.setRe(re);
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

#ifdef OUTPUT
	geom_gen.print();
#endif

	if (geom_type == 0){
		std::string filenameParam("../uq_data/uq_parameter_" + std::to_string(id) + ".param");
		geom_gen.setFilenames("../uq_data/uq_geometry.geom", filenameParam.c_str());
	}

	geom_gen.writeToFile();
}

int main(int argc, char** argv){
	//frontend_commandlineargs(argc,argv);
	frontend_interactive();
	return 0;
}

