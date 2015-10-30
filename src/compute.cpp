#include "compute.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "solver.hpp"

#include <cmath>	// sin, M_PI
#include <algorithm>    // std::min
#include <iostream>	// std::cout

/* Public methods */

/* Constructor */
Compute::Compute(const Geometry *geom, const Parameter *param)
: _t(0.0), _dtlimit(0.0), _epslimit(0.0)
{
	// TODO: Werte fuer _dtlimit, _epslimit richtig?
	//_epslimit = param->Eps();
	_epslimit = 1e-4;
	//_dtlimit = 

	_geom = geom;
	_param = param;
	//_u = new Grid(_geom, multi_real_t(_geom->Mesh()[0],_geom->Mesh()[1]/2.0)); // offset here?
	//_v = new Grid(_geom, multi_real_t(_geom->Mesh()[0]/2.0,_geom->Mesh()[1]));
	//_p = new Grid(_geom, multi_real_t(_geom->Mesh()[0]/2.0,_geom->Mesh()[1]/2.0));
	_u = new Grid(_geom, multi_real_t(-1.0,-0.5));
	_v = new Grid(_geom, multi_real_t(-0.5,-1.0));
	_p = new Grid(_geom, multi_real_t(-0.5,-0.5));
	_F = new Grid(_geom);
	_G = new Grid(_geom);
	_rhs = new Grid(_geom);
	//_tmp = new Grid(_geom);

	/*Grid u(_geom);
	u.Initialize(0.0);*/ // just for debugging issues

	// initialize grids
	//std::cout << "Compute: Initializing the grids..." << std::flush; // only for debugging issues
	//std::cout << "\n_u at adress " << _u  << "\n" << std::flush; // only for debugging issues
	_u->Initialize(0.0);
	_v->Initialize(0.0);

	_p->Initialize(0.0);

	_F->Initialize(0.0);
	_G->Initialize(0.0);
	_rhs->Initialize(0.0);
	//_temp->Initialize(0.0);

	//set boundary values
	update_boundary_values();

	/*Iterator it(_geom);
	it.First();
	while (it.Valid()){
		if (it.Pos()[1]>=64){
			_p->Cell(it) = 10.01;
		}
		//std::cout << it.Pos()[0] << ", " << it.Pos()[1] << "\n";
		it.Next();
	}*/

	//std::cout << "Done.\n" << std::flush; // only for debugging issues

	// write boundary values
	//std::cout << "Compute: Updating the boundary values..." << std::flush; // only for debugging issues
	//update_boundary_values();
	//std::cout << "Done.\n" << std::flush; // only for debugging issues
	
	// TODO: is this right, anything else to initialize?
	real_t h = 0.5 * (_geom->Mesh()[0] + _geom->Mesh()[1]); // just took the average here
	real_t omega = 2.0 / (1.0+sin(M_PI*h)); // TODO: set omega to the right value
	//real_t omega = 1.0;
	
	_solver = new SOR(_geom, omega);
	//_solver = new JacobiSolver(_geom);

	//_solver = new HeatConductionSolver(_geom);
}

/* Destructor */
Compute::~Compute()
{
	// TODO: something more to delete?
	delete[] _u;
	delete[] _v;
	delete[] _p;
	delete[] _F;
	delete[] _G;
	delete[] _rhs;
	delete[] _tmp;
	
	delete _solver;
}

void Compute::TimeStep(bool printInfo, bool verbose=false)
{
	// TODO: test
	
	// compute dt
	//if (verbose) std::cout << "Computing the timestep width..." << std::flush; // only for debugging issues
	real_t dt = compute_dt();
	//if (verbose) std::cout << "Done.\n" << std::flush; // only for debugging issues
	
	// boundary values
	// update_boundary_values();
	
	// compute F, G
	MomentumEqu(dt);
	
	// compute rhs
	RHS(dt);
	//_rhs->Initialize(12.0);
	//_solver->delete_average(_rhs);

	//_geom->Update_V(_p); // TODO: remove
	
	// solve Poisson equation
	real_t residual(_epslimit + 1.0);
	index_t iteration(0);
	//while (iteration <= _param->IterMax() && residual > _epslimit){
	while (true){
		// do one solver cycle here
		
		// boundary values
		//update_boundary_values();
		
		//_geom->Update_P(_p);

		residual = _solver->Cycle(_p, _rhs);

		// delete average
		//_solver->delete_average(_p);

		iteration++;
		if (iteration > _param->IterMax()){
			std::cout << "Warning: Solver did not converge! Residual: " << residual << "\n";
			break;
		} else if (residual < _epslimit){
			std::cout << "Solver converged after " << iteration << " timesteps. Residual: " << residual << "\n";
			break;
		}
	}
	
	// compute u, v
	NewVelocities(dt);
	update_boundary_values();
	_u->CheckNaN();
	_v->CheckNaN();

	//update total time
	_t += dt;

	// print information
	if (printInfo){
		std::cout << "============================================================\n";
		// total simulated time
		std::cout << "Total simulated time: t = " << _t << "\n";
		// timestep
		std::cout << "Last timestep: dt = " << dt << "\n";
		// magnitudes of the fields
		std::cout << "max(F) = " << _F->AbsMax() << ", max(G) = " << _G->AbsMax() << ", max(rhs) = " << _rhs->AbsMax() << "\n";
		std::cout << "max(u) = " << _u->AbsMax() << ", max(v) = " << _v->AbsMax() << ", max(p) = " << _p->AbsMax() << "\n";
		std::cout << "Average value of rhs: " << _rhs->average_value() << "\n";
	}
}

const real_t& Compute::GetTime() const
{
	return _t;
}

const Grid* Compute::GetU() const
{
	return _u;
}

const Grid* Compute::GetV() const
{
	return _v;
}

const Grid* Compute::GetP() const
{
	return _p;
}

const Grid* Compute::GetRHS() const
{
	return _rhs;
}

const Grid* Compute::GetVelocity()
{
	// TODO: test
	//Grid* res = new Grid(_geom,multi_real_t((-1)*_geom->Mesh()[0],(-1)*_geom->Mesh()[1]));
	Grid* res = new Grid(_geom,multi_real_t(-0.5,-0.5));	
	Iterator it(_geom);
	it.First();
	while (it.Valid()){
		res->Cell(it) = sqrt(pow(0.5*(_u->Cell(it)+_u->Cell(it.Left())),2.0)+pow(0.5*(_v->Cell(it)+_v->Cell(it.Down())),2.0));
		it.Next();
	}
	res->CheckNaN();
	return res;
}

const Grid* Compute::GetVorticity()
{
	// TODO: test
	Grid* res = new Grid(_geom);
	InteriorIterator it(_geom);
	it.First();
	while (it.Valid()){
		// here, we use central difference quotients
		res->Cell(it) = _v->dx_c(it) - _u->dy_c(it);
		it.Next();
	}
	return res;
}

const Grid* Compute::GetStream()
{
	// TODO
}

/* private methods */

void Compute::NewVelocities(const real_t& dt)
{
	// TODO: test
	InteriorIterator it(_geom);
	it.First();
	while (it.Valid()){
		_u->Cell(it) = _F->Cell(it) - dt * _p->dx_r(it);
		_v->Cell(it) = _G->Cell(it) - dt * _p->dy_r(it);
		it.Next();
	}
}

void Compute::MomentumEqu(const real_t& dt)
{
	// TODO: test
	InteriorIterator it(_geom);
	it.First();
	while (it.Valid()){
		_F->Cell(it) = _u->Cell(it) + dt * ( 1.0/_param->Re() * (_u->dxx(it) + _u->dyy(it)) - _u->DC_duu_x(it,_param->Alpha()) - _u->DC_duv_y(it,_param->Alpha(),_v) );
		_G->Cell(it) = _v->Cell(it) + dt * ( 1.0/_param->Re() * (_v->dxx(it) + _v->dyy(it)) - _v->DC_dvv_y(it,_param->Alpha()) - _v->DC_duv_x(it,_param->Alpha(),_u) );
		it.Next();
	}
}

void Compute::RHS(const real_t& dt)
{
	// TODO: test
	InteriorIterator it(_geom);
	it.First();
	while (it.Valid()){
		_rhs->Cell(it) = 1.0/dt * (_F->dx_l(it) + _G->dy_l(it));
		it.Next();
	}
	//_solver->delete_average(_rhs);
}

// own methods
real_t Compute::compute_dt() const
{
	//std::cout << "max u: " << _u->AbsMax() << ", max v: " << _v->AbsMax() << "\n"; // for debugging issues
	real_t res = std::min(_geom->Mesh()[0] / _u->AbsMax(), _geom->Mesh()[1] / _v->AbsMax());
	res = std::min(res, _param->Re()/2.0 * pow(_geom->Mesh()[0],2.0) * pow(_geom->Mesh()[1],2.0) / (pow(_geom->Mesh()[0],2.0)+pow(_geom->Mesh()[1],2.0)));
	res *= _param->Tau(); // just to be sure (because it is a strict inequality)
	//res /= 25.0;
	return res;
}

void Compute::update_boundary_values()
{
	_geom->Update_U(_u);
	_geom->Update_V(_v);

	_geom->Update_U(_F);
	_geom->Update_V(_G);

	//_geom->Update_P(_p); // is this the correct place for this update?
}
