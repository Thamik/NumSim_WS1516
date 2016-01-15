#include "parameter.hpp"
#include "typedef.hpp"
#include "communicator.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>     /* atof */

#include <mpi.h>

/**
The constructor of the Parameter class. Constructs a default parameter set.
*/
Parameter::Parameter(Communicator* comm)
: _re(1000.0), _omega(1.7), _alpha(0.9), _eps(0.001), _tau(0.5), _itermax(100), _itermin(0), _gravity(0.0,-9.81), _pr(3.0), _beta(0.5), _gamma(0.9), _q(1.0), _comm(comm) // standard parameter
{
	// First time parameter values (see exercise 1, page 7 below)
	_dt = 0.2;
	_tend = 82 * 0.2;
}

/**
This method loads a parameter set from a file. In this file, the parameters should be listed in different rows in the following order: Reynolds number, relaxation factor, upwind differencing factor, tolerance for pressure iteration, safety factor for the time step size, maximum number of iterations, time step size, end time.
\param[in]	file	Filename of the file, which contains the parameter set
\param[in]	verbose	Determines, if further information shall be printed to the standard output channel (usually the command line)
*/
void Parameter::Load(const char* file, bool verbose)
{
for (int i_rank=0; i_rank < _comm->getSize(); i_rank++){
if(i_rank == _comm->getRank()){ // read the files sequentially

	if (verbose && _comm->getRank()==0){
		std::cout << "Loading parameter file from path " << file << " ...\n";
	}
	std::string temp_string;
	std::ifstream infile;
	infile.open(file);
	if (!infile.is_open()){
		std::cout << "Warning: parameter file could not be read!\n";
		continue;
	}
	for(int i=1;i<=9;i++)
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
				_itermin = atoi(temp_string.c_str());
				break;
			case 8:
				_dt = atof(temp_string.c_str());
				break;
			case 9:
				_tend = atof(temp_string.c_str());
				break;
		}
	}
	infile.close();
	//if (verbose){
#ifdef OUTPUT_PARAMS
	if (_comm->getRank() == 0){
		std::cout << "--------------------------------------------------\n";
		std::cout << "Parameter configuration read from file:\n";
		std::cout << "Re\t\t=\t" << _re << "\n";
		std::cout << "Omega\t\t=\t" << _omega << "\n";
		std::cout << "Alpha\t\t=\t" << _alpha << "\n";
		std::cout << "Epsilon\t\t=\t" << _eps << "\n";
		std::cout << "Tau\t\t=\t" << _tau << "\n";
		std::cout << "Itermax\t\t=\t" << _itermax << "\n";
		std::cout << "Itermin\t\t=\t" << _itermin << "\n";
		std::cout << "dt\t\t=\t" << _dt << "\n";
		std::cout << "tend\t\t=\t" << _tend << "\n";
		std::cout << "Gravity\t\t=\t(" << _gravity[0] << ", " << _gravity[1] << ")\n";
		std::cout << "Pr\t\t=\t" << _pr << "\n";
		std::cout << "Beta\t\t=\t" << _beta << "\n";
		std::cout << "Gamma\t\t=\t" << _gamma << "\n";
		std::cout << "Q\t\t=\t" << _q << "\n";
		std::cout << "--------------------------------------------------\n";
	}
#endif

}
MPI_Barrier(MPI_COMM_WORLD);
}

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

const index_t& Parameter::IterMin() const
{
	return _itermin;
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

const multi_real_t& Parameter::Gravity() const
{
	return _gravity;
}

const real_t& Parameter::Pr() const
{
	return _pr;
}

const real_t& Parameter::Beta() const
{
	return _beta;
}

const real_t& Parameter::Gamma() const
{
	return _gamma;
}

const real_t& Parameter::Q() const
{
	return _q;
}

