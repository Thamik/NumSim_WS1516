#include "parameter.hpp"
#include "typedef.hpp"

/* Constructor */
Parameter::Parameter()
: _re(1000.0), _omega(1.7), _alpha(0.9), _eps(0.001), _tau(0.5), _itermax(100)
{
	// TODO: richtige Werte fuer dt und tend
	_dt = 0.1; // was muss hier eigentlich rein?
	_tend = 10.0;
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
