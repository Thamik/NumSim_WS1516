#include "solver.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "grid.hpp"

#include <cmath>	// std::abs
#include <algorithm>	// std::max

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
	// TODO: test
	return std::abs(rhs->Cell(it) - grid->dxx(it) - grid->dyy(it));
}

// own method
real_t Solver::totalRes(const Grid* grid, const Grid* rhs) const
{
	// TODO: test
	real_t totalRes(0.0);
	InteriorIterator it(_geom);
	it.First();
	while (it.Valid()){
		totalRes = std::max(totalRes, localRes(it,grid,rhs));
		it.Next();
	}
	return totalRes;
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
	InteriorIterator it(_geom);
	it.First();
	while (it.Valid()){
		//corr = pow(_geom->Mesh()[0],2.0) * pow(_geom->Mesh()[1],2.0) / (2.0 * (pow(_geom->Mesh()[0],2.0)+pow(_geom->Mesh()[1],2.0))) * ( grid->dxx(it) + grid->dyy(it) - (2.0 * (pow(_geom->Mesh()[0],2.0)+pow(_geom->Mesh()[1],2.0))) * grid->Cell(it) / (pow(_geom->Mesh()[0],2.0) * pow(_geom->Mesh()[1],2.0)) - rhs->Cell(it) );
		//correct correction term (note, that dxx(it) and dyy(it) calculate wrong fractions (they include the middle term!))
		corr = pow(_geom->Mesh()[0],2.0) * pow(_geom->Mesh()[1],2.0) / (2.0 * (pow(_geom->Mesh()[0],2.0)+pow(_geom->Mesh()[1],2.0))) * ( grid->dxx(it) + grid->dyy(it) - rhs->Cell(it) );
		//corr = localRes(it, grid, rhs);

		grid->Cell(it) += _omega * corr;
		//grid->Cell(it) = (1.0-_omega) * grid->Cell(it) + _omega * corr;

		it.Next();
	}
	return totalRes(grid,rhs);
}
