#include "parameter.hpp"
#include "typedef.hpp"

/* Constructor */
Parameter::Parameter()
{
	// TODO: hier fehlen noch die default Parameter
}

void Parameter::Load(const char* file)
{
	// TODO: hier fehlt noch alles
}

/* Getter functions */
const real_t& Parameter::Re() const
{
	return _re;
}

const real_t& Parameter::Omega() const
{
	return _omega;
}

const real_t& Parameter::Alpha() const
{
	return _alpha;
}

const real_t& Parameter::Dt() const
{
	return _dt;
}

const real_t& Parameter::Tend() const
{
	return _tend;
}

const index_t& Parameter::IterMax() const
{
	return _itermax;
}

const real_t& Parameter::Eps() const
{
	return _eps;
}

const real_t& Parameter::Tau() const
{
	return _tau;
}
