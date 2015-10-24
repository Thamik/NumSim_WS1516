#include "compute.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include "geometry.hpp"
#include "parameter.hpp"

#include <cmath>
#include <algorithm>    // std::min

/* Public methods */

/* Constructor */
Compute::Compute(const Geometry *geom, const Parameter *param)
: _t(0.0), _dtlimit(0.0), _epslimit(0.0)
{
	// TODO: Werte fuer _dtlimit, _epslimit richtig?

	_geom = geom;
	_param = param;
	_u = new Grid(_geom); // offset here?
	_v = new Grid(_geom);
	_p = new Grid(_geom);
	_F = new Grid(_geom);
	_G = new Grid(_geom);
	_rhs = new Grid(_geom);
	_tmp = new Grid(_geom);
	
	// TODO: initialize solver etc.
}

/* Destructor */
Compute::~Compute()
{
	// TODO: something more to delete?
	delete _u;
	delete _v;
	delete _p;
	delete _F;
	delete _G;
	delete _rhs;
	delete _tmp;
}

void Compute::TimeStep(bool printInfo)
{
	// TODO
	
	// compute dt
	real_t dt = compute_dt();
	
	// boundary values
	
	// compute F, G
	MomentumEqu(dt);
	
	// compute rhs
	
	// solve Poisson equation
	
	// compute u, v
	
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
	Grid* res = new Grid(_geom);
	Iterator it(_geom);
	it.First();
	while (it.Valid()){
		res->Cell(it) = sqrt(pow(_u->Cell(it),2.0)+pow(_v->Cell(it),2.0));
		it.Next();
	}
	return res;
}

const Grid* Compute::GetVorticity()
{
	// TODO: test
	Grid* res = new Grid(_geom);
	Iterator it(_geom);
	it.First();
	while (it.Valid()){
		// here, we use central difference quotients
		res->Cell(it) = 0.5*(_v->dx_l(it)+_v->dx_r(it)) - 0.5*(_u->dy_l(it)+_u->dy_r(it));
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
	// TODO
}

void Compute::MomentumEqu(const real_t& dt)
{
	// TODO
}

void Compute::RHS(const real_t& dt)
{
	// TODO
}

// own methods
real_t Compute::compute_dt() const
{
	real_t res = std::min(_geom->Mesh()[0] / _u->AbsMax(), _geom->Mesh()[1] / _v->AbsMax());
	res = std::min(res, _param->Re()/2.0 * pow(_geom->Mesh()[0],2.0) * pow(_geom->Mesh()[1],2.0) / (pow(_geom->Mesh()[0],2.0)+pow(_geom->Mesh()[1],2.0)));
	res /= 2.0; // just to be sure (because it is a strict inequality)
	return res;
}
