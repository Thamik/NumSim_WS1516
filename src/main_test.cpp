#include "typedef.hpp"
#include "compute.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "solver.hpp"
#include "visu.hpp"
#include "vtk.hpp"

#include <iostream>

int main(int argc, char **argv)
{
	Parameter param;
	Geometry geom;
	SOR solver(&geom,1.0);
	//JacobiSolver solver(&geom);
	
	Grid p(&geom);
	Grid laplace(&geom);
	laplace.Initialize(0.0);
	p.Initialize(0.0);

	Grid rhs(&geom);
	rhs.Initialize(2.0);
	rhs.Out();

	Renderer visu(geom.Length(), geom.Mesh());
	visu.Init(800,800);
	visu.Render(&p);

	p.Out();
	std::cout << "Size of geom:" << geom.Size()[0] << "x" << geom.Size()[1];
	std::cin.get();
	
	real_t res(0.0);
	for (int i=1; i<=100; i++){

		// TODO: right order of: solve - BC Update - delete_av.??
		res = solver.Cycle(&p, &rhs);
		//Neumann RB
		geom.Update_P(&p);
		//delete mean value
		solver.delete_average(&p);


		visu.Render(&p);
		std::cout << "Residual: " << res << ", max: " << p.AbsMax() << "\n";

		// Output of p
		p.Out();
		
		// Calculating and printing of laplacian of p
		laplace.Laplace(&p);
		geom.Update_P(&laplace);		
		laplace.Out();
		std::cin.get(); // pause
	}
	

	return 0;
}
