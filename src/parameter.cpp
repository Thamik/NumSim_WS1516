#include "parameter.hpp"

/* Constructor */
Parameter::Parameter()
{
	// hier fehlen noch die default Parameter
}

void Parameter::Load(const char *file)
{
	// hier fehlt noch alles
}

/* Getter functions */
const real_t& Parameter::Re()
{
	return _re;
}

const real_t& Parameter::Omega()
{
	return _omega;
}

const real_t& Parameter::Alpha()
{
	return _alpha;
}

const real_t& Parameter::Dt()
{
	return _dt;
}

const real_t& Parameter::Tend()
{
	return _tend;
}

const index_t& Parameter::IterMax()
{
	return _itermax;
}

const real_t& Parameter::Eps()
{
	return _eps;
}

const real_t& Parameter::Tau()
{
	return _tau;
}
