#include "compute.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "solver.hpp"
#include "communicator.hpp"

#include <cmath>	// sin, M_PI
#include <algorithm>    // std::min
#include <iostream>	// std::cout

/* Public methods */

/* Constructor */
/**
\param[in]	geom	pointer on geometry object providing the geometry information
\param[in]	param	pointer on parameter object providing the parameter data
\param[in]	comm 	pointer on Communicator object
*/
Compute::Compute(const Geometry *geom, const Parameter *param, const Communicator *comm)
: _t(0.0), _dtlimit(0.0), _epslimit(0.0), _geom(geom), _param(param), _comm(comm)
{
	// TODO: Werte fuer _dtlimit, _epslimit richtig?
	_epslimit = param->Eps();
	//_epslimit = 1e-4;
	//_dtlimit = 

	_u = new Grid(_geom, multi_real_t(-1.0,-0.5));
	_v = new Grid(_geom, multi_real_t(-0.5,-1.0));
	_p = new Grid(_geom, multi_real_t(-0.5,-0.5));
	_F = new Grid(_geom);
	_G = new Grid(_geom);
	_rhs = new Grid(_geom);

	_tmp_velocity = new Grid(_geom, multi_real_t(-0.5,-0.5));
	_tmp_vorticity = new Grid(_geom);
	_tmp_stream = new Grid(_geom);

	// Initialize Grids
	_u->Initialize(0.0);
	_v->Initialize(0.0);

	_p->Initialize(0.0);

	_F->Initialize(0.0);
	_G->Initialize(0.0);
	_rhs->Initialize(0.0);
	
	_tmp_velocity->Initialize(0.0);
	_tmp_vorticity->Initialize(0.0);
	_tmp_stream->Initialize(0.0);

	//set boundary values
	update_boundary_values();

	// TODO: is this right, anything else to initialize?
	real_t h = 0.5 * (_geom->Mesh()[0] + _geom->Mesh()[1]); // just took the average here
	real_t omega = 2.0 / (1.0+sin(M_PI*h)); // TODO: set omega to the right value
	//real_t omega = 1.0;
	
	_solver = new RedOrBlackSOR(_geom, omega);
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

	delete _tmp_velocity;
	delete _tmp_vorticity;
	delete _tmp_stream;
	
	delete _solver;
}

/**
In order to do so, first a reasonable dt for stability is calculated and F, G, RHS are evaluated.
Afterwards, the Poisson pressure equation is solved and the velocitys are updated.
\param[in] printInfo boolean if additional informations on the fields and rediduum of p are printed
\param[in] verbose boolean if debbuging information should be printed (standard: false)
*/
void Compute::TimeStep(bool printInfo, bool verbose)
{
	// TODO: test
	// compute dt
	if (verbose) std::cout << "Computing the timestep width..." << std::flush; // only for debugging issues
	real_t dt = compute_dt(); // BLOCKING
	if (verbose) std::cout << "Done.\n" << std::flush; // only for debugging issues
	
	// compute F, G...
	MomentumEqu(dt);
	update_boundary_values(); //update boundary values

	// ...and sync them
	sync_FG(); // BLOCKING

	// compute rhs
	RHS(dt);
	
	// solve Poisson equation
	real_t residual(_epslimit + 1.0);
	index_t iteration(0);
	while (true){
		// one solver cycle is done here, alternating red or black
		if ((iteration % 2) == 0){
			residual = _solver->RedCycle(_p, _rhs);
		} else {
			residual = _solver->BlackCycle(_p, _rhs);
		}

		sync_p();

		iteration++;

		if ((iteration/2) > _param->IterMax()){ // iteration/2 because we only do half a cycle in each iteration
			if (_comm->getRank()==0) std::cout << "Warning: Solver did not converge! Residual: " << residual << "\n";
			break;
		} else if (residual < _epslimit){
			if (_comm->getRank()==0) std::cout << "Solver converged after " << iteration << " iterations. Residual: " << residual << "\n";
			break;
		}
	}
	
	// compute new velocitys u, v...
	NewVelocities(dt);
	update_boundary_values();

	// ...and sync them
	sync_uv(); // BLOCKING

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
		//std::cout << "max(F) = " << _F->TotalAbsMax() << ", max(G) = " << _G->TotalAbsMax() << ", max(rhs) = " << _rhs->TotalAbsMax() << "\n";
		//std::cout << "max(u) = " << _u->TotalAbsMax() << ", max(v) = " << _v->TotalAbsMax() << ", max(p) = " << _p->TotalAbsMax() << "\n";
		
		//std::cout << "Average value of rhs: " << _rhs->average_value() << "\n";
		std::cout << "============================================================\n";
	}

	// print debug information
	/*if (_comm->getRank() == 3) std::cout << "(" << _comm->getRank() << "a): " << _u->Data()[16+18*2] << ", " << _v->Data()[16+18*2] << "\n" << std::flush;
	if (_comm->getRank() == 3) std::cout << "(" << _comm->getRank() << "b): " << _u->Data()[16+18] << ", " << _v->Data()[16+18] << "\n" << std::flush;
	if (_comm->getRank() == 1) std::cout << "(" << _comm->getRank() << "a): " << _u->Data()[18*17-2] << ", " << _v->Data()[18*17-2] << "\n" << std::flush;
	if (_comm->getRank() == 1) std::cout << "(" << _comm->getRank() << "b): " << _u->Data()[18*16-2] << ", " << _v->Data()[18*16-2] << "\n" << std::flush;*/
}

/**
\return actual total time
*/
const real_t& Compute::GetTime() const
{
	return _t;
}

/**
\return pointer on the Grid containing u
*/
const Grid* Compute::GetU() const
{
	return _u;
}

/**
\return pointer on the Grid containing v
*/
const Grid* Compute::GetV() const
{
	return _v;
}

/**
\return pointer on the Grid containing p
*/
const Grid* Compute::GetP() const
{
	return _p;
}

/**
\return pointer on the Grid containing the RHS
*/
const Grid* Compute::GetRHS() const
{
	return _rhs;
}

/**
The absolute velocity in the euklidean norm is calculted at the middle point of the physical grid cells
\return pointer on the Grid containing absolute velocity
*/
const Grid* Compute::GetVelocity()
{
	// TODO: test
	Iterator it(_geom);
	it.First();
	while (it.Valid()){
		_tmp_velocity->Cell(it) = sqrt(pow(0.5*(_u->Cell(it)+_u->Cell(it.Left())),2.0)+pow(0.5*(_v->Cell(it)+_v->Cell(it.Down())),2.0));
		it.Next();
	}
	//_tmp_velocity->CheckNaN();
	return _tmp_velocity;
}

/**
\return pointer on the Grid containing the Vorticity
*/
const Grid* Compute::GetVorticity()
{
	// TODO: test
	//InteriorIterator it(_geom);
	Iterator it(_geom);
	it.First();
	while (it.Valid()){
		// here, we use central difference quotients
		//_tmp_vorticity->Cell(it) = _v->dx_c(it) - _u->dy_c(it);
		
		// as wanted in the exercise sheet, we use forward quotients
		_tmp_vorticity->Cell(it) = _u->dy_r(it) - _v->dx_r(it);

		it.Next();
	}
	return _tmp_vorticity;
}

/**
\return pointer on the Grid containing the stream
*/
const Grid* Compute::GetStream()
{
	// set values in own grid
	Iterator it(_geom);
	it.First();
	_tmp_stream->Cell(it) = 0;
	it.Next();
	while (it.Valid()){
		if (it.Pos()[1] == 0){
			_tmp_stream->Cell(it) = _tmp_stream->Cell(it.Left()) - _v->Cell(it) * _geom->Mesh()[0];
		} else {
			_tmp_stream->Cell(it) = _tmp_stream->Cell(it.Down()) + _u->Cell(it) * _geom->Mesh()[1];
		}
		it.Next();
	}

	// communicate and update values at the bottom subdomains
	index_t nx = _comm->ThreadDim()[0];
	index_t ny = _comm->ThreadDim()[1];
	for (int i=0; i<=(int(nx)-2); i++){
		if (_comm->getRank() == _comm->getRankDistribution(i,0)){
			MPI_Send(&(_tmp_stream->Data()[_geom->Size()[0]-1-2]), 1, MPI_DOUBLE, _comm->getRankDistribution(i+1, 0), 0, MPI_COMM_WORLD);
		} else if (_comm->getRank() == _comm->getRankDistribution(i+1, 0)){
			real_t recBuf(0);
			MPI_Recv(&recBuf, 1, MPI_DOUBLE, _comm->getRankDistribution(i, 0), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			// update values (add constant)
			Iterator it(_geom);
			it.First();
			_tmp_stream->Cell(it) = recBuf - _v->Cell(it) * _geom->Mesh()[0];
			it.Next();
			while (it.Valid()){
				_tmp_stream->Cell(it) += _tmp_stream->Data()[0];
				it.Next();
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}

	// comunicate and update all other values
	for (int j=0; j<=(int(ny)-2); j++){
		for (int i=0; i<int(nx); i++){
			if (_comm->getRank() == _comm->getRankDistribution(i, j)){
				MPI_Send(&(_tmp_stream->Data()[_geom->Size()[0]*(_geom->Size()[1]-1-2)]), 1, MPI_DOUBLE, _comm->getRankDistribution(i, j+1), 0, MPI_COMM_WORLD);
			} else if (_comm->getRank() == _comm->getRankDistribution(i, j+1)){
				real_t recBuf(0);
				MPI_Recv(&recBuf, 1, MPI_DOUBLE, _comm->getRankDistribution(i, j), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// update values (add constant)
				Iterator it(_geom);
				it.First();
				_tmp_stream->Cell(it) = recBuf + _u->Cell(it) * _geom->Mesh()[1];
				it.Next();
				while (it.Valid()){
					_tmp_stream->Cell(it) += _tmp_stream->Data()[0];
					it.Next();
				}
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}

	//if (_comm->getRank()==0) std::cout << "Stream has been computed!\n" << std::flush;

	return _tmp_stream;
}

/* private methods */

/**
The new velocitys are calculated and stored in the _u,_v grids (in-place
\param[in] dt time step size
*/
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

/**
The expression for F,G arrising for the Momentum equations is evaluated and stored in the Grids _F,_G
\param[in] dt time step size
*/
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

/**
The Right-Hand-Side of the pressure poisson equation is calculated depending on F,G and stored into _rhs
\param[in] dt time step size
*/
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
/**
In the algorithm, dt has to be sufficiently small to provide stability of the algorithm. This is garanteed by calculating the minimum and multiplying with a safety factor
\return stabile time step
*/
real_t Compute::compute_dt() const
{
	real_t u_absmax = _u->TotalAbsMax();
	real_t v_absmax = _v->TotalAbsMax();
	real_t res = std::min(_geom->Mesh()[0]/u_absmax, _geom->Mesh()[1]/v_absmax);
	res = std::min(res, _param->Re()/2.0 * pow(_geom->Mesh()[0],2.0) * pow(_geom->Mesh()[1],2.0) / (pow(_geom->Mesh()[0],2.0)+pow(_geom->Mesh()[1],2.0)));
	res *= _param->Tau(); // just to be sure (because it is a strict inequality)
	return res;
}

void Compute::update_boundary_values()
{
	_geom->Update_U(_u);
	_geom->Update_V(_v);

	_geom->Update_U(_F);
	_geom->Update_V(_G);
}

void Compute::sync_FG()
{
	_comm->copyBoundary(_F);
	_comm->copyBoundary(_G);
}

void Compute::sync_uv()
{
	_comm->copyBoundary(_u);
	_comm->copyBoundary(_v);
}

void Compute::sync_p()
{
	_comm->copyBoundary(_p);
}

void Compute::sync_all()
{
	sync_FG();
	sync_uv();
	sync_p();
}
