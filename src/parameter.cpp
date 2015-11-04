#include "parameter.hpp"
#include "typedef.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>     /* atof */

/**
The constructor of the Parameter class. Constructs a default parameter set.
*/
Parameter::Parameter()
: _re(1000.0), _omega(1.7), _alpha(0.9), _eps(0.001), _tau(0.5), _itermax(100) // standard parameter
{
	// First time parameter values (see exercise 1, page 7 below)
	_dt = 0.2;
	_tend = 82 * 0.2;
}

/**
This method loads a parameter set from a file. In this file, the parameters should be listed in different rows in the following order: Reynolds number, relaxation factor, upwind differencing factor, tolerance for pressure iteration, safety factor for the time step size, maximum number of iterations, time step size, end time.
\param[in]	file	Filename of the file, which contains the parameter set
*/
void Parameter::Load(const char* file)
{
	// TODO: Test this method
	std::string temp_string;
	std::ifstream infile;
	infile.open(file);
	for(int i=1;i<=8;i++)
	{
		getline(infile,temp_string);
		switch(i)
		{
			case 1:
				_re = atof(temp_string.c_str());
				break;
			case 2:
				_omega = atof(temp_string.c_str());
				break;
			case 3:
				_alpha = atof(temp_string.c_str());
				break;
			case 4:
				_eps = atof(temp_string.c_str());
				break;
			case 5:
				_tau = atof(temp_string.c_str());
				break;
			case 6:
				_itermax = atoi(temp_string.c_str());
				break;
			case 7:
				_dt = atof(temp_string.c_str());
				break;
			case 8:
				_tend = atof(temp_string.c_str());
				break;
		}
	}
	infile.close();
}

/**
Getter function for the Reynolds number.
\return	The Reynolds number
*/
const real_t& Parameter::Re() const
{
	return _re;
}

/**
Getter function for the relaxation factor.
\return	The relaxation factor omega
*/
const real_t& Parameter::Omega() const
{
	return _omega;
}

/**
Getter function for the upwind differencing factor.
\return	The upwind differencing factor alpha
*/
const real_t& Parameter::Alpha() const
{
	return _alpha;
}

/**
Getter function for the time step size.
\return	The time step size
*/
const real_t& Parameter::Dt() const
{
	return _dt;
}

/**
Getter function for the end time.
\return	The end time
*/
const real_t& Parameter::Tend() const
{
	return _tend;
}

/**
Getter function for the maximum number of iterations.
\return	The maximum number of iterations
*/
const index_t& Parameter::IterMax() const
{
	return _itermax;
}

/**
Getter function for the tolerance for the pressure iteration.
\return	The tolerance epsilon
*/
const real_t& Parameter::Eps() const
{
	return _eps;
}

/**
Getter function for the safety factor for the time step size.
\return	The safety factor tau
*/
const real_t& Parameter::Tau() const
{
	return _tau;
}
