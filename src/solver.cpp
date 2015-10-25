#include "solver.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "grid.hpp"

#include <cmath>

// Solver

// Constructor
Solver::Solver(const Geometry* geom)
: _geom(geom)
{
}

// Destructor
Solver::~Solver()
{
}

real_t Solver::localRes(const Iterator& it, const Grid* grid, const Grid* rhs) const
{
	// TODO
}

// Concrete SOR solver

// Constructor
SOR::SOR(const Geometry* geom, const real_t& omega)
: Solver(geom), _omega(omega)
{
}

// Destructor
SOR::~SOR()
{
}

real_t SOR::Cycle(Grid* grid, const Grid* rhs) const
{
	// TODO: test
	real_t corr(0.0);
	Iterator it(_geom);
	it.First();
	while (it.Valid()){
		corr = pow(_geom->Mesh()[0],2.0) * pow(_geom->Mesh()[1],2.0) / (2.0 * (pow(_geom->Mesh()[0],2.0)+pow(_geom->Mesh()[1],2.0))) * ( grid->dxx(it) + grid->dyy(it) - (2.0 * (pow(_geom->Mesh()[0],2.0)+pow(_geom->Mesh()[1],2.0))) * grid->Cell(it) / (pow(_geom->Mesh()[0],2.0) * pow(_geom->Mesh()[1],2.0)) - rhs->Cell(it) );
		grid->Cell(it) += _omega * corr;
		it.Next();
	}
}
