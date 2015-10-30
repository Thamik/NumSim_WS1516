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

	//std::cout << "rhs in localRes: " << rhs->Cell(it) << "\n";
	
	//return std::abs(- grid->dxx(it) - grid->dyy(it));
}

// own method
real_t Solver::totalRes(const Grid* grid, const Grid* rhs) const
{
	// TODO: test
	real_t totalRes(0.0);
	InteriorIterator it(_geom);
	it.First();
	while (it.Valid()){
		//std::cout << "localRes " << it.Pos()[0] << ", " << it.Pos()[1] << ": " << localRes(it,grid,rhs) << "\n";
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
	//grid->Cell(it.Down()) = 0.0;
	while (it.Valid()){
		//correct correction term (note, that dxx(it) and dyy(it) calculate wrong fractions (they include the middle term!))
		real_t hx = _geom->Mesh()[0] * _geom->Mesh()[0];
		real_t hy = _geom->Mesh()[1] * _geom->Mesh()[1];

		//corr = pow(_geom->Mesh()[0],2.0) * pow(_geom->Mesh()[1],2.0) / (2.0 * (pow(_geom->Mesh()[0],2.0)+pow(_geom->Mesh()[1],2.0))) * ( grid->dxx(it) + grid->dyy(it) - rhs->Cell(it) );

		corr = (hx*hy)/(2.0*(hx + hy)) * ( grid->dxx(it) + grid->dyy(it) - rhs->Cell(it) );
		grid->Cell(it) += _omega * corr;

		//Neuman BC TODO
		if (it.Right().Right().Value() == it.Right().Value())
			grid->Cell(it.Right()) = grid->Cell(it);
		if (it.Left().Left().Value() == it.Left().Value())
			grid->Cell(it.Left()) = grid->Cell(it);
		if (it.Down().Down().Value() == it.Down().Value())
			grid->Cell(it.Down()) = grid->Cell(it);
		if (it.Top().Top().Value() == it.Top().Value())
			grid->Cell(it.Top()) = grid->Cell(it);

		
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
		
		grid->Cell(it) = 1.0/(-2.0/pow(_geom->Mesh()[0],2.0) - 2.0/pow(_geom->Mesh()[1],2.0)) * (rhs->Cell(it) - 1.0/pow(_geom->Mesh()[0],2.0) * cpy->Cell(it.Left()) - 1.0/pow(_geom->Mesh()[0],2.0) * cpy->Cell(it.Right()) - 1.0/pow(_geom->Mesh()[1],2.0) * cpy->Cell(it.Top()) - 1.0/pow(_geom->Mesh()[0],2.0) * cpy->Cell(it.Down()));

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


HeatConductionSolver::HeatConductionSolver(const Geometry* geom)
: Solver(geom)
{
}

HeatConductionSolver::~HeatConductionSolver()
{
}

real_t HeatConductionSolver::Cycle(Grid* grid, const Grid* rhs) const
{
	real_t dt = 0.000001;
	Grid* cpy = grid->copy();
	for (int i=0; i<=10; i++){
		InteriorIterator it(_geom);
		it.First();
		while (it.Valid()){
			grid->Cell(it) += dt * (cpy->dxx(it)+cpy->dyy(it) - rhs->Cell(it));
			it.Next();
		}
	}
	delete cpy;
	return totalRes(grid,rhs);
}

