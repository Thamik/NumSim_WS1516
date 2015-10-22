#include "compute.hpp"

/* Public methods */

/* Constructor */
Compute::Compute(const Geometry *geom, const Parameter *param)
{
	// TODO
}

/* Destructor */
Compute::~Compute()
{
	// TODO
}

void Compute::TimeStep(bool printInfo)
{
	// TODO
}

const real_t& Compute::GetTime() const
{
	// TODO
}

const Grid* Compute::GetU() const
{
	// TODO
}

const Grid* Compute::GetV() const
{
	// TODO
}

const Grid* Compute::GetP() const
{
	// TODO
}

const Grid* Compute::GetRHS() const
{
	// TODO
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
