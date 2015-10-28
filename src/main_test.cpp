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
	JacobiSolver solver(&geom);
	
	Grid p(&geom);
	p.Initialize(0.0);

	Grid rhs(&geom);
	rhs.Initialize(2.0);

	Renderer visu(geom.Length(), geom.Mesh());
	visu.Init(800,800);
	visu.Render(&p);
	std::cin.get();
	
	real_t res(0.0);
	for (int i=1; i<=100; i++){
		geom.Update_P(&p);
		res = solver.Cycle(&p, &rhs);
		visu.Render(&p);
		std::cout << "Residual: " << res << ", max: " << p.AbsMax() << "\n";
		std::cin.get(); // pause
	}
	

	return 0;
}
