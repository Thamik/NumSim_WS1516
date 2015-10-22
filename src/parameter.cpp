#include "parameter.hpp"
#include "typedef.hpp"

//#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>     /* atof */

/* Constructor */
Parameter::Parameter()
: _re(1000.0), _omega(1.7), _alpha(0.9), _eps(0.001), _tau(0.5), _itermax(100)
{
	// Erste Zeit-Parameter-Werte (siehe Blatt 1, Seite 7 unten)
	_dt = 0.2;
	_tend = 82 * 0.2;
}

void Parameter::Load(const char* file)
{
	// TODO: Test this method
	std::string temp_string;
	ifstream infile;
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
