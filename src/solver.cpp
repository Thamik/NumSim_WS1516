#include "solver.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "grid.hpp"

#include <cmath>	// std::abs
#include <algorithm>	// std::max
#include <iostream>

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
		//if (localRes(it,grid,rhs) < 1.99) std::cout << it.Pos()[0] << ", " << it.Pos()[1] << "\n";
	}
	return totalRes;
}

void Solver::delete_average(Grid* grid) const
{
	// compute the average value
	real_t avg(0.0);
	Iterator it(_geom);
	it.First();
	while (it.Valid()){
		avg += grid->Cell(it);
		it.Next();
	}
	avg /= real_t(_geom->Size()[0] * _geom->Size()[1]);
	Iterator it2(_geom);
	it2.First();
	while (it2.Valid()){
		grid->Cell(it2) -= avg;
		//std::cout << "1: " << grid->Data()[it2.Value()] << "\n";
		//grid->Data()[it2.Value()] -= avg;
		//std::cout << "2: " << grid->Data()[it2.Value()] << "\n";
		it2.Next();
	}
	//std::cout << "Average: " << avg << "\n";
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


JacobiSolver::JacobiSolver(const Geometry* geom)
: Solver(geom)
{
}

JacobiSolver::~JacobiSolver()
{
}

real_t JacobiSolver::Cycle(Grid* grid, const Grid* rhs) const
{
	Grid* cpy = grid->copy();
	InteriorIterator it(_geom);
	it.First();
	while (it.Valid()){
		grid->Cell(it) = 1.0/(-2.0/_geom->Mesh()[0] - 2.0/_geom->Mesh()[1]) * (rhs->Cell(it) - 1.0/_geom->Mesh()[0] * cpy->Cell(it.Left()) - 1.0/_geom->Mesh()[0] * cpy->Cell(it.Right()) - 1.0/_geom->Mesh()[1] * cpy->Cell(it.Top()) - 1.0/_geom->Mesh()[0] * cpy->Cell(it.Down()));

		/*if (std::abs(grid->Cell(it) - cpy->Cell(it))>1e-8){
			std::cout << "hier hat sich was geändert: " << it.Pos()[0] << ", " << it.Pos()[1] << "\n";
		} else {
			std::cout << "diff: " << 1.0/(-2.0/_geom->Mesh()[0] - 2.0/_geom->Mesh()[1]) * (rhs->Cell(it) - 1.0/_geom->Mesh()[0] * cpy->Cell(it.Left()) - 1.0/_geom->Mesh()[0] * cpy->Cell(it.Right()) - 1.0/_geom->Mesh()[1] * cpy->Cell(it.Top()) - 1.0/_geom->Mesh()[0] * cpy->Cell(it.Down())) - grid->Cell(it) << "\n";
		}*/

		it.Next();
	}
	delete cpy;
	return totalRes(grid,rhs);
}
