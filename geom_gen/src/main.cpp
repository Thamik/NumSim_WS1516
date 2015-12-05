#include <iostream>

#include "geom_gen.hpp"

int main(int argc, char** argv){

	GeometryGenerator geom_gen;

	geom_gen.setSize(60,10);
	//geom_gen.drivenCavity();
	//geom_gen.pipeFlow();
	geom_gen.flowOverAStep();

	geom_gen.writeToFile("../data/complex_default.geom");

}
