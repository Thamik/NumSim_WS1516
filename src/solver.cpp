#include "solver.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "grid.hpp"
#include "communicator.hpp"

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
#ifdef COMPLEX_GEOM
	InteriorIteratorGG it(_geom);
#else
	InteriorIterator it(_geom);
#endif
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
#ifdef COMPLEX_GEOM
	InteriorIteratorGG it(_geom);
#else
	InteriorIterator it(_geom);
#endif
	it.First();
	while (it.Valid()){
		totalRes += localRes(it,grid,rhs);
		it.Next();
	}
	return totalRes/((_geom->Size()[0]-1.0) * (_geom->Size()[1]-1.0));
}

/**
\param[in]	grid	Approximate solution
\param[in]	rhs	Right hand side of the pressure Poisson equation
\return		The (synchronized) total residual of all grids
*/
real_t Solver::synced_totalRes(const Grid* grid, const Grid* rhs) const
{
	//return synced_totalRes_L1_averaged(grid,rhs);
	return synced_totalRes_Linf(grid,rhs);
}

/**
\param[in]	grid	Approximate solution
\param[in]	rhs	Right hand side of the pressure Poisson equation
\return		The (synchronized) total Linf residual of all grids
*/
real_t Solver::synced_totalRes_Linf(const Grid* grid, const Grid* rhs) const
{
	return _geom->getCommunicator()->gatherMax(totalRes_Linf(grid,rhs));
}

/**
\param[in]	grid	Approximate solution
\param[in]	rhs	Right hand side of the pressure Poisson equation
\return		The (synchronized) total L1 residual of all grids
*/
real_t Solver::synced_totalRes_L1_averaged(const Grid* grid, const Grid* rhs) const
{
	return (_geom->getCommunicator()->gatherSum(totalRes_L1_averaged(grid,rhs))) / (_geom->getCommunicator()->getSize());
}

/**
This method deletes the average value from an arbitrary scalar grid.
\param[in][out]	grid	A scalar grid, which average value is being deleted
*/
void Solver::delete_average(Grid* grid) const
{
	// compute the average value
	real_t avg = grid->average_value();

#ifdef COMPLEX_GEOM
	InteriorIteratorGG it(_geom);
#else
	InteriorIterator it(_geom);
#endif

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
#ifdef COMPLEX_GEOM
	InteriorIteratorGG it(_geom);
#else
	InteriorIterator it(_geom);
#endif
	it.First();
	while (it.Valid()){
		erase_local_residual(grid,rhs,it);
		it.Next();
	}
	// update the pressure boundary values
	_geom->Update_P(grid);
	return totalRes(grid,rhs);
}

/**
\param[in][out]	grid	Approximate solution
\param[in]	rhs	Right hand side of the pressure Poisson equation
\param[in]	it	Iterator which specifies the current position in the grid
*/
void SOR::erase_local_residual(Grid* grid, const Grid* rhs, const Iterator& it) const
{
	real_t hx = _geom->Mesh()[0] * _geom->Mesh()[0];
	real_t hy = _geom->Mesh()[1] * _geom->Mesh()[1];

	real_t corr = (hx*hy)/(2.0*(hx + hy)) * ( grid->dxx(it) + grid->dyy(it) - rhs->Cell(it) );
	grid->Cell(it) += _omega * corr;

	/*if (grid->Cell(it) < -0.2) {
		std::cout << it.Pos()[0] << ", " << it.Pos()[1] << std::endl;
		std::cout << rhs->Cell(it) << std::endl;
	}*/

	// set one cell to zero (the alternative to deleting the pressure average)
	//if (it.Value()==((_geom->Size()[0]*_geom->Size()[1])/2)+_geom->Size()[0]/2){
	/*if (it.Value()==((_geom->TotalSize()[0]*_geom->TotalSize()[1])/2)+_geom->TotalSize()[0]/2){
		grid->Cell(it) = 0.0;
	}*/
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
This is the constructor of the RedOrBlackSOR class.
\param[in] 	geom 	A Geometry object that encapsulates all necessary geometrical information
\param[in]	omega	Relaxation factor
*/
RedOrBlackSOR::RedOrBlackSOR(const Geometry* geom, const real_t& omega)
: SOR(geom, omega)
{
	// TODO: something else to do here?
}

/**
This is the destructor of the RedOrBlackSOR class.
*/
RedOrBlackSOR::~RedOrBlackSOR()
{
	// TODO: is there something to do here?
}

/**
Performs a red cycle to solve the Poisson equation.
\param[in][out]	grid	Initial value and result of the cycle
\param[in]	rhs	Grid, where the right hand side of the Poisson equation is specified
*/
real_t RedOrBlackSOR::RedCycle(Grid* grid, const Grid* rhs) const
{
#ifdef COMPLEX_GEOM
	JumpingInteriorIteratorGG it(_geom, false);
#else
	JumpingInteriorIterator it(_geom, false);
#endif
	it.First();

	real_t hx = _geom->Mesh()[0] * _geom->Mesh()[0];
	real_t hy = _geom->Mesh()[1] * _geom->Mesh()[1];
	real_t twohxhyinv = 1.0/(2.0*(hx+hy));

	while (it.Valid()){
//		erase_local_residual(grid,rhs,it);

		real_t corr = ( (grid->Cell(it.Left())-2.0*grid->Cell(it)+grid->Cell(it.Right()))*hy + (grid->Cell(it.Top())-2.0*grid->Cell(it)+grid->Cell(it.Down()))*hx - rhs->Cell(it)*hx*hy) * twohxhyinv;
		grid->Cell(it) += _omega * corr;

		it.Next();
	}
	// update the pressure boundary values
#ifdef COMPLEX_GEOM
	_geom->UpdateGG_P(grid);
#else
	_geom->Update_P(grid);
#endif
	return synced_totalRes(grid,rhs);
}

/**
Performs a black cycle to solve the Poisson equation.
\param[in][out]	grid	Initial value and result of the cycle
\param[in]	rhs	Grid, where the right hand side of the Poisson equation is specified
*/
real_t RedOrBlackSOR::BlackCycle(Grid* grid, const Grid* rhs) const
{
#ifdef COMPLEX_GEOM
	JumpingInteriorIteratorGG it(_geom, true);
#else
	JumpingInteriorIterator it(_geom, true);
#endif
	it.First();

	real_t hx = _geom->Mesh()[0] * _geom->Mesh()[0];
	real_t hy = _geom->Mesh()[1] * _geom->Mesh()[1];
	real_t twohxhyinv = 1.0/(2.0*(hx+hy));

	while (it.Valid()){
//		erase_local_residual(grid,rhs,it);

		real_t corr = ( (grid->Cell(it.Left())-2.0*grid->Cell(it)+grid->Cell(it.Right()))*hy + (grid->Cell(it.Top())-2.0*grid->Cell(it)+grid->Cell(it.Down()))*hx - rhs->Cell(it)*hx*hy) * twohxhyinv;
		grid->Cell(it) += _omega * corr;

		it.Next();
	}
	// update the pressure boundary values
#ifdef COMPLEX_GEOM
	_geom->UpdateGG_P(grid);
#else
	_geom->Update_P(grid);
#endif
	return synced_totalRes(grid,rhs);
}

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
/*real_t JacobiSolver::Cycle(Grid* grid, const Grid* rhs) const
{
	// TODO: implement the update of the boundary pressure data
	Grid* cpy = grid->copy();
	InteriorIteratorGG it(_geom);
	it.First();
	while (it.Valid()){
		grid->Cell(it) = 1.0/(-2.0/pow(_geom->Mesh()[0],2.0) - 2.0/pow(_geom->Mesh()[1],2.0)) * (rhs->Cell(it) - 1.0/pow(_geom->Mesh()[0],2.0) * cpy->Cell(it.Left()) - 1.0/pow(_geom->Mesh()[0],2.0) * cpy->Cell(it.Right()) - 1.0/pow(_geom->Mesh()[1],2.0) * cpy->Cell(it.Top()) - 1.0/pow(_geom->Mesh()[0],2.0) * cpy->Cell(it.Down()));

		it.Next();
	}
	delete cpy;
	return totalRes(grid,rhs);
}*/
// the parallel version
real_t JacobiSolver::Cycle(Grid* grid, const Grid* rhs) const
{
	Grid* cpy = grid->copy();
#ifdef COMPLEX_GEOM
	InteriorIteratorGG it(_geom);
#else
	InteriorIterator it(_geom);
#endif
	it.First();
	while (it.Valid()){
		grid->Cell(it) = 1.0/(-2.0/pow(_geom->Mesh()[0],2.0) - 2.0/pow(_geom->Mesh()[1],2.0)) * (rhs->Cell(it) - 1.0/pow(_geom->Mesh()[0],2.0) * cpy->Cell(it.Left()) - 1.0/pow(_geom->Mesh()[0],2.0) * cpy->Cell(it.Right()) - 1.0/pow(_geom->Mesh()[1],2.0) * cpy->Cell(it.Top()) - 1.0/pow(_geom->Mesh()[0],2.0) * cpy->Cell(it.Down()));

		it.Next();
	}
	delete cpy;
#ifdef COMPLEX_GEOM
	_geom->UpdateGG_P(grid);
#else
	_geom->Update_P(grid);
#endif
	return synced_totalRes(grid,rhs);
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
#ifdef COMPLEX_GEOM
		InteriorIteratorGG it(_geom);
#else
		InteriorIterator it(_geom);
#endif
		it.First();
		while (it.Valid()){
			grid->Cell(it) += dt * (cpy->dxx(it)+cpy->dyy(it) - rhs->Cell(it));
			it.Next();
		}
	}
	delete cpy;
	return totalRes(grid,rhs);
}

//------------------------------------------------------------------------------

/*MGInfoHandle::MGInfoHandle()
{
}

MGInfoHandle::~MGInfoHandle()
{
}*/

MGSolver::MGSolver(real_t eps, index_t itermax)
: _eps(eps), _itermax(itermax)
{
}

MGSolver::~MGSolver()
{
}

real_t MGSolver::Solve(Grid* pressure, const Grid* rhs) const
{
	real_t total_res(_eps * 2.0);
	bool converged(false);

	// compute coarser geometry
	Geometry geom_coarse(pressure->getGeometry()->getCommunicator());
	geom_coarse.halfSize(pressure->getGeometry());

	index_t minGridSize = 16;

	if (std::min(geom_coarse.TotalSize()[0]-2, geom_coarse.TotalSize()[1]-2) <= minGridSize){
		// do not use a coarser grid anymore, solve per SOR

		real_t h = 0.5 * (pressure->getGeometry()->Mesh()[0] + pressure->getGeometry()->Mesh()[1]);
		real_t omega = 2.0 / (1.0+sin(M_PI*h));
		RedOrBlackSOR solver = RedOrBlackSOR(pressure->getGeometry(), omega);
		index_t iteration(0);
		while(true){
			if ((iteration % 2) == (pressure->getGeometry()->getCommunicator()->EvenOdd() ? 1 : 0)){
				total_res = solver.RedCycle(pressure, rhs);
			} else {
				total_res = solver.BlackCycle(pressure, rhs);
			}

			if (pressure->getGeometry()->getCommunicator()->getSize() > 1){
				pressure->getGeometry()->getCommunicator()->copyBoundary(pressure);
			}

			iteration++;

			if (iteration > _itermax*2){
				converged = false;
				break;
			} else if (total_res < _eps){
				converged = true;
				break;
			}
		}
	} else {
		// solve on coarser grid

		// allocate memory for the grids
		Grid res(pressure->getGeometry(),multi_real_t(-0.5,-0.5));
		Grid res_coarse(&geom_coarse,multi_real_t(-0.5,-0.5));
		Grid err_coarse(&geom_coarse,multi_real_t(-0.5,-0.5));
		Grid err(pressure->getGeometry(),multi_real_t(-0.5,-0.5));

		// do W-cycle
		index_t iter(0);
		index_t max_w_cycles = 50;
		while (total_res >= _eps && iter <= max_w_cycles){

			// smooth
			total_res = smooth(pressure, rhs);

			// compute residual
			//res.Initialize(0.0); // this should not be needed
			compute_residual(pressure, rhs, &res);

			// restrict to coarser grid
			//res_coarse.Initialize(0.0); // this should not be needed
			restrict_grid(&res, &res_coarse);

			// MG iteration
			err_coarse.Initialize(0.0); // THIS IS NEEDED
			Solve(&err_coarse, &res_coarse); // solve recursively with homogeneous boundary

			// prolongate to finer grid
			err.Initialize(0.0); // this should not be needed, but seems to be
			interpolate_grid(&err_coarse, &err);

			// ... and finally add the error to this grid
			pressure->add(&err);

			// postsmooting
			total_res = smooth(pressure, rhs);

			iter++;

		}
	}
	
	return total_res;
}

real_t MGSolver::smooth(Grid* pressure, const Grid* rhs) const
{
	real_t res(0.0);

	RedOrBlackSOR solver = RedOrBlackSOR(pressure->getGeometry(), 1.0); // use Gauss-Seidel solver
	index_t n_cycles = 5;
	for (index_t ii=0; ii<n_cycles*2; ii++){
		if ((ii % 2) == (pressure->getGeometry()->getCommunicator()->EvenOdd() ? 1 : 0)){
			res = solver.RedCycle(pressure, rhs);
		} else {
			res = solver.BlackCycle(pressure, rhs);
		}

		if (pressure->getGeometry()->getCommunicator()->getSize() > 1){
			pressure->getGeometry()->getCommunicator()->copyBoundary(pressure);
		}
	}

	return res;
}

void MGSolver::compute_residual(const Grid* pressure, const Grid* rhs, Grid* res) const
{
	//InteriorIteratorGG it(pressure->getGeometry());
	Iterator it(pressure->getGeometry());
	it.First();
	while(it.Valid()) {
		res->Cell(it) = rhs->Cell(it) - pressure->dxx(it) - pressure->dyy(it);
		it.Next();
	}
}

void MGSolver::restrict_grid(const Grid* old_grid, Grid* new_grid) const
{
	InteriorIteratorGG it_new(new_grid->getGeometry());
	it_new.First();

	/*while(it_new.Valid()) {
		// multi_real_t physpos = it_new.PhysPos(new_grid->getOffset());
		multi_real_t physpos = it_new.PhysPos(multi_real_t(-0.5,-0.5));
		if (physpos[0] < 0 || physpos[0] > 1 || physpos[1] < 0 || physpos[1] > 1){
			std::cout << "Warning! Physpos: " << physpos[0] << ", " << physpos[1] << std::endl;
		}
		new_grid->Cell(it_new) = old_grid->Interpolate(physpos);
		it_new.Next();
	}*/

	multi_index_t size = new_grid->getGeometry()->Size();
	index_t posX, posY, x1, x2, y1, y2, vallu, valru, vallo, valro, newsizeX;

	newsizeX = 2*size[0]-2;

	while(it_new.Valid()) {
		posX = it_new.Pos()[0];
		posY = it_new.Pos()[1];

		x1 = 2*posX-1;
		x2 = 2*posX;
		y1 = 2*posY-1;
		y2 = 2*posY;

		//calculate values
		vallu = x1 + y1*newsizeX;
		valru = x2 + y1*newsizeX;
		vallo = x1 + y2*newsizeX;	
		valro = x2 + y2*newsizeX;

		//interpolation
		new_grid->Cell(it_new) = 0.25 * (old_grid->Data()[vallu] + old_grid->Data()[valru] + old_grid->Data()[vallo] + old_grid->Data()[valro]);

		it_new.Next();
	}

	new_grid->getGeometry()->UpdateGG_P(new_grid);
}

void MGSolver::interpolate_grid(const Grid* old_grid, Grid* new_grid) const
{
	/*InteriorIteratorGG it_new(new_grid->getGeometry());
	it_new.First();

	while(it_new.Valid()) {
		//multi_real_t physpos = it_new.PhysPos(new_grid->getOffset());
		multi_real_t physpos = it_new.PhysPos(multi_real_t(-0.5,-0.5));
		if (physpos[0] < 0 || physpos[0] > 1 || physpos[1] < 0 || physpos[1] > 1){
			std::cout << "Warning! Physpos: " << physpos[0] << ", " << physpos[1] << std::endl;
		}
		new_grid->Cell(it_new) = old_grid->Interpolate(physpos);
		it_new.Next();
	}*/

	InteriorIteratorGG it_old(old_grid->getGeometry());
	it_old.First();

	multi_index_t size = old_grid->getGeometry()->Size();
	index_t posX, posY, x1, x2, y1, y2, vallu, valru, vallo, valro, newsizeX;

	newsizeX = 2*size[0]-2;

	while(it_old.Valid()) {
		posX = it_old.Pos()[0];
		posY = it_old.Pos()[1];

		x1 = 2*posX-1;
		x2 = 2*posX;
		y1 = 2*posY-1;
		y2 = 2*posY;

		//calculate values
		vallu = x1 + y1*newsizeX;
		valru = x2 + y1*newsizeX;
		vallo = x1 + y2*newsizeX;	
		valro = x2 + y2*newsizeX;

		new_grid->Data()[vallu] = old_grid->Cell(it_old);
		new_grid->Data()[valru] = old_grid->Cell(it_old);
		new_grid->Data()[vallo] = old_grid->Cell(it_old);
		new_grid->Data()[valro] = old_grid->Cell(it_old);

		it_old.Next();
	}

	new_grid->getGeometry()->UpdateGG_P(new_grid);
}

