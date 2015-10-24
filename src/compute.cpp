#include "compute.hpp"
#include "grid.hpp"

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
	// TODO
}

const Grid* Compute::GetVorticity()
{
	// TODO
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
