#include "solver.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "grid.hpp"

#include <cmath>	// std::abs
#include <algorithm>	// std::max
#include <iostream>

//------------------------------------------------------------------------------

/* Solver */

/**
The constructor of the abstract Solver class.
\param[in]	geom	Geometry object providing all necessary geometrical information
*/
Solver::Solver(const Geometry* geom)
: _geom(geom)
{
}

/**
The destructor of the abstract Solver class.
*/
Solver::~Solver()
{
}

/**
This method returns the local residual for the pressure Poisson equation specified by grid and rhs at the position given by it.
\param[in]	it	Iterator, which specifies the position, where the local residual shall be computed
\param[in]	grid	Approximate solution of the pressure Poisson equation
\param[in]	rhs	Right hand side of the pressure Poisson equation

\return		The local residual at [it]
*/
real_t Solver::localRes(const Iterator& it, const Grid* grid, const Grid* rhs) const
{
	return std::abs(rhs->Cell(it) - grid->dxx(it) - grid->dyy(it));
}

/*real_t Solver::localRes(const Iterator& it, const Grid* grid, const Grid* rhs) const
{
	index_t eps_W(0);
	index_t eps_E(0);
	index_t eps_S(0);
	index_t eps_N(0);

	if (it.Pos()[0]==1) eps_W = 0; else eps_W = 1;
	if (it.Pos()[0]==_geom->Size()[0]-2) eps_E = 0; else eps_E = 1;
	if (it.Pos()[1]==1) eps_S = 0; else eps_S = 1;
	if (it.Pos()[1]==_geom->Size()[1]-2) eps_N = 0; else eps_N = 1;

	real_t hx = _geom->Mesh()[0] * _geom->Mesh()[0];
	real_t hy = _geom->Mesh()[1] * _geom->Mesh()[1];

	return pow((eps_E*(grid->Cell(it.Right())-grid->Cell(it)) - eps_W*(grid->Cell(it)-grid->Cell(it.Left())))/hx + (eps_N*(grid->Cell(it.Top())-grid->Cell(it)) - eps_S*(grid->Cell(it)-grid->Cell(it.Down())))/hy - rhs->Cell(it),2.0);
}*/

/**
This function return the total residual of the approximate solution in grid for the discrete pressure Poisson equation with right hand side like specified in rhs.
\param[in]	grid	Approximate solution
\param[in]	rhs	Right hand side of the pressure Poisson equation

\return 	The total residual
*/
real_t Solver::totalRes(const Grid* grid, const Grid* rhs) const
{
	return totalRes_L1_averaged(grid,rhs);
}

/**
This function return the total Linf residual of the approximate solution in grid for the discrete pressure Poisson equation with right hand side like specified in rhs.
\param[in]	grid	Approximate solution
\param[in]	rhs	Right hand side of the pressure Poisson equation

\return 	The total Linf residual
*/
real_t Solver::totalRes_Linf(const Grid* grid, const Grid* rhs) const
{
	real_t totalRes(0.0);
	InteriorIterator it(_geom);
	it.First();
	while (it.Valid()){
		totalRes = std::max(totalRes, localRes(it,grid,rhs));
		it.Next();
	}
	return totalRes;
}

/**
This function return the averaged L1 residual of the approximate solution in grid for the discrete pressure Poisson equation with right hand side like specified in rhs.
\param[in]	grid	Approximate solution
\param[in]	rhs	Right hand side of the pressure Poisson equation

\return 	The averaged L1 residual
*/
real_t Solver::totalRes_L1_averaged(const Grid* grid, const Grid* rhs) const
{
	real_t totalRes(0.0);
	InteriorIterator it(_geom);
	it.First();
	while (it.Valid()){
		totalRes += localRes(it,grid,rhs);
		it.Next();
	}
	return totalRes/((_geom->Size()[0]-1.0) * (_geom->Size()[1]-1.0));
}

/**
This method deletes the average value from an arbitrary scalar grid.
\param[in][out]	grid	A scalar grid, which average value is being deleted
*/
void Solver::delete_average(Grid* grid) const
{
	// compute the average value
	real_t avg = grid->average_value();

	Iterator it(_geom);
	it.First();
	while (it.Valid()){
		grid->Cell(it) -= avg;
		it.Next();
	}
}

//------------------------------------------------------------------------------

/* Concrete SOR solver */

/**
The constructor of the SOR class.
\param[in]	geom	Geoemetry object containing all geometrical information needed
\param[in]	omega	Relaxation factor
*/
SOR::SOR(const Geometry* geom, const real_t& omega)
: Solver(geom), _omega(omega)
{
}

/**
The destructor of the SOR class.
*/
SOR::~SOR()
{
}

/**
This method executes a SOR solver cycle on the given grid and returns the total residual.
\param[in][out]	grid	The current pressure values
\param[in]	rhs	The right hand side of the pressure Poisson equation

\return 	The total residual of the approximate solution in grid after the solver cycle
*/
real_t SOR::Cycle(Grid* grid, const Grid* rhs) const
{
	real_t corr(0.0);
	InteriorIterator it(_geom);
	it.First();
	while (it.Valid()){
		real_t hx = _geom->Mesh()[0] * _geom->Mesh()[0];
		real_t hy = _geom->Mesh()[1] * _geom->Mesh()[1];

		corr = (hx*hy)/(2.0*(hx + hy)) * ( grid->dxx(it) + grid->dyy(it) - rhs->Cell(it) );
		grid->Cell(it) += _omega * corr;

		// set one cell to zero (the alternative to deleting the pressure average)
		if (it.Value()==((_geom->Size()[0]*_geom->Size()[1])/2)+_geom->Size()[0]/2){
			grid->Cell(it) = 0.0;
		}

		it.Next();
	}
	// update the pressure boundary values
	_geom->Update_P(grid);
	return totalRes(grid,rhs);
}

// SOR Cycle, Neumann boundary conditions included
// see: Numerical Simulation in Fluid Dynamics, p.37
/*real_t SOR::Cycle(Grid* grid, const Grid* rhs) const
{
	// TODO: test

	index_t eps_W(0);
	index_t eps_E(0);
	index_t eps_S(0);
	index_t eps_N(0);

	InteriorIterator it(_geom);
	it.First();
	//grid->Cell(it.Down()) = 0.0;
	while (it.Valid()) {

		if (it.Pos()[0]==1) eps_W = 0; else eps_W = 1;
		if (it.Pos()[0]==_geom->Size()[0]-2) eps_E = 0; else eps_E = 1;
		if (it.Pos()[1]==1) eps_S = 0; else eps_S = 1;
		if (it.Pos()[1]==_geom->Size()[1]-2) eps_N = 0; else eps_N = 1;

		real_t hx = _geom->Mesh()[0] * _geom->Mesh()[0];
		real_t hy = _geom->Mesh()[1] * _geom->Mesh()[1];

		grid->Cell(it) = (1-_omega)*grid->Cell(it) + _omega / ( (eps_E+eps_W)/hx + (eps_N+eps_S)/hy ) * ( (eps_E * grid->Cell(it.Right()) + eps_W * grid->Cell(it.Left()))/hx + (eps_N * grid->Cell(it.Top()) + eps_S * grid->Cell(it.Down()))/hy - rhs->Cell(it) );

		// set one cell to zero
		if (it.Value()==((_geom->Size()[0]*_geom->Size()[1])/2)+_geom->Size()[0]/2){
			grid->Cell(it) = 0.0;
		}
		
		it.Next();
	}

	_geom->Update_P(grid);

	return totalRes(grid,rhs);
}*/

//------------------------------------------------------------------------------

/* Red Or Black SOR Solver */

/**
This is the constructor of the RedOrBlackSOR Solver class.
\param[in] 	geom 	A Geometry object that 

//------------------------------------------------------------------------------

/* Jacobi Solver */

/**
This is the constructor of the JacobiSolver class.
\param[in]	geom	A Geometry object that encapsulates all necessary geometrical information
*/
JacobiSolver::JacobiSolver(const Geometry* geom)
: Solver(geom)
{
}

/**
This is the destructor of the JacobiSolver class.
*/
JacobiSolver::~JacobiSolver()
{
}

/**
This method executes a Jacobi solver cycle on the given grid and returns the total residual.
\param[in][out]	grid	The current pressure values
\param[in]	rhs	The right hand side of the pressure Poisson equation

\return 	The total residual of the approximate solution in grid after the solver cycle
*/
real_t JacobiSolver::Cycle(Grid* grid, const Grid* rhs) const
{
	// TODO: implement the update of the boundary pressure data
	Grid* cpy = grid->copy();
	InteriorIterator it(_geom);
	it.First();
	while (it.Valid()){
		grid->Cell(it) = 1.0/(-2.0/pow(_geom->Mesh()[0],2.0) - 2.0/pow(_geom->Mesh()[1],2.0)) * (rhs->Cell(it) - 1.0/pow(_geom->Mesh()[0],2.0) * cpy->Cell(it.Left()) - 1.0/pow(_geom->Mesh()[0],2.0) * cpy->Cell(it.Right()) - 1.0/pow(_geom->Mesh()[1],2.0) * cpy->Cell(it.Top()) - 1.0/pow(_geom->Mesh()[0],2.0) * cpy->Cell(it.Down()));

		it.Next();
	}
	delete cpy;
	return totalRes(grid,rhs);
}

//------------------------------------------------------------------------------

/* Heat Conduction Solver */

/**
This is the constructor of the HeatConductionSolver class.
\param[in]	geom	A Geometry object that specifies all geometrical data needed
*/
HeatConductionSolver::HeatConductionSolver(const Geometry* geom)
: Solver(geom)
{
}

/**
This is the destructor method of the HeatConductionSolver class.
*/
HeatConductionSolver::~HeatConductionSolver()
{
}

/**
This method executes a solver cycle on the given grid and returns the total residual. To do this, the solution of the corresponding heat conduction equation is approximated and the limit of this solution for the time converging towards infinity is computed.
\param[in][out]	grid	The current pressure values
\param[in]	rhs	The right hand side of the pressure Poisson equation

\return 	The total residual of the approximate solution in grid after the solver cycle
*/
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
