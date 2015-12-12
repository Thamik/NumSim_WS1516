#include "compute.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "solver.hpp"
#include "communicator.hpp"
#include "console_output.hpp"
#include "particles.hpp"

#include <cmath>	// sin, M_PI
#include <algorithm>    // std::min
#include <iostream>	// std::cout

#include <math.h>

/* Public methods */

/* Constructor */
/**
\param[in]	geom	pointer on geometry object providing the geometry information
\param[in]	param	pointer on parameter object providing the parameter data
\param[in]	comm 	pointer on Communicator object
*/
Compute::Compute(const Geometry *geom, const Parameter *param, const Communicator *comm)
: _t(0.0), _dtlimit(0.0), _epslimit(0.0), _solver_converging(false), _diff_p(0.0), _diff_rhs(0.0), _diff_F(0.0), _diff_G(0.0), _geom(geom), _param(param), _comm(comm), _particles(nullptr)
{
	// TODO: Werte fuer _dtlimit, _epslimit richtig?
	_epslimit = param->Eps();
	//_epslimit = 1e-4;
	//_dtlimit = 

	if (_comm->getRank() == 0){
		// build particles object
		_particles = new Particles(_geom);

		// define options
//		_particles->streaklinePolicy();
		_particles->particleTracingPolicy();
//		_particles->setMatlabFormat();
//		_particles->setMatlabOneFileFormat();
		_particles->setPythonOneFileFormat();

		_particles->init();
	}

	// construct grids
	_u = new Grid(_geom, multi_real_t(-1.0,-0.5));
	_v = new Grid(_geom, multi_real_t(-0.5,-1.0));
	_p = new Grid(_geom, multi_real_t(-0.5,-0.5));
	_F = new Grid(_geom);
	_G = new Grid(_geom);
	_rhs = new Grid(_geom);
	_tmp_velocity = new Grid(_geom, multi_real_t(-0.5,-0.5));
	_tmp_vorticity = new Grid(_geom);
	_tmp_stream = new Grid(_geom);

	_p_old = new Grid(_geom, multi_real_t(-0.5,-0.5));
	_rhs_old = new Grid(_geom);
	_F_old = new Grid(_geom);
	_G_old = new Grid(_geom);

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

	_p_old->Initialize(0.0);
	_rhs_old->Initialize(0.0);
	_F_old->Initialize(0.0);
	_G_old->Initialize(0.0);

	//set boundary values
	update_boundary_values();

	// construct solver
	real_t h = 0.5 * (_geom->Mesh()[0] + _geom->Mesh()[1]); // just took the average here
	real_t omega = 2.0 / (1.0+sin(M_PI*h));
	//real_t omega = 1.0;	
	_solver = new RedOrBlackSOR(_geom, omega);
	//_solver = new JacobiSolver(_geom);

	// construct console clock
	_clock = new ConsoleClock();
}

/* Destructor */

/**
The destructor of the Compute class.
*/
Compute::~Compute()
{
	delete _u;
	delete _v;
	delete _p;
	delete _F;
	delete _G;
	delete _rhs;

	delete _tmp_velocity;
	delete _tmp_vorticity;
	delete _tmp_stream;

	delete _p_old;
	delete _rhs_old;
	delete _F_old;
	delete _G_old;
	
	delete _solver;

	delete _clock;

	if (_comm->getRank() == 0 && _particles != nullptr){
		_particles->finalizeFile();
		delete _particles;
	}
}

/**
In order to do so, first a reasonable dt for stability is calculated and F, G, RHS are evaluated.
Afterwards, the Poisson pressure equation is solved and the velocitys are updated.
\param[in] printInfo boolean if additional informations on the fields and rediduum of p are printed
\param[in] verbose boolean if debbuging information should be printed (standard: false)
*/
void Compute::TimeStep(bool verbose, real_t diff_time)
{
	// compute dt
	if (verbose) std::cout << "Computing the timestep width..." << std::flush; // only for debugging issues
	real_t dt = compute_dt(); // BLOCKING

	bool printInfo(false);
	if (diff_time <= dt){
		dt = diff_time;
		printInfo = !_comm->getRank();
	}
	printInfo = !_comm->getRank();
	//printInfo = false;

	/*if (dt < 0){
		std::cout << "Fatal error! Compute::TimeStep(): negative timestep!" << std::endl;
		throw std::runtime_error("compute: negative timestep");
	} else if (dt < 1e-7){
		std::cout << "Warning! Compute::TimeStep(): very small timestep!" << std::endl;
	} else if (dt > 1){
		std::cout << "Warning! Compute::TimeStep(): very large timestep!" << std::endl;
	}*/

	if (verbose) std::cout << "Done.\n" << std::flush; // only for debugging issues
	
	// for debugging issues
	delete _F_old;
	delete _G_old;
	_F_old = _F->copy();
	_G_old = _G->copy();

	// compute F, G...
	MomentumEqu(dt);
	update_boundary_values(); //update boundary values

	// ...and sync them
	sync_FG(); // BLOCKING

	// for debugging issues
	delete _rhs_old;
	_rhs_old = _rhs->copy();

	// compute rhs
	RHS(dt);

	//_geom->UpdateGG_P(_p);

	//update_boundary_values(); // debug

	// for debugging issues
	delete _p_old;
	_p_old = _p->copy();
	
	// solve Poisson equation
	real_t residual(_epslimit + 1.0);
	index_t iteration(0);
	while (true){
		// one solver cycle is done here, alternating red or black
		if ((iteration % 2) == (_comm->EvenOdd() ? 1 : 0)){
			residual = _solver->RedCycle(_p, _rhs);
		} else {
			residual = _solver->BlackCycle(_p, _rhs);
		}

		// for Jacobi solvers
		//residual = _solver->Cycle(_p, _rhs);

		sync_p();

		iteration++;

		if (iteration >= _param->IterMin()){
			// check for edge condition
			if ((iteration/2) > _param->IterMax()){ // iteration/2 because we only do half a cycle in each iteration
				//if (_comm->getRank()==0) std::cout << "Warning: Solver did not converge! Residual: " << residual << "\n";
				_solver_converging = false;
				break;
			} else if (residual < _epslimit){
				//if (_comm->getRank()==0) std::cout << "Solver converged after " << iteration << " iterations. Residual: " << residual << "\n";
				_solver_converging = true;
				break;
			}
		}
	}
	_geom->UpdateGG_P(_p);

	update_boundary_values(); // debug

	check_for_incontinuities();
	
	// compute new velocitys u, v...
	NewVelocities(dt);
	update_boundary_values();

	// ...and sync them
	sync_uv(); // BLOCKING

	update_boundary_values(); // debug

	//update total time
	_t += dt;

	//update infos
	currentIterations += iteration;
	currentTime += dt;
	currentNoTimeSteps++;
	currentResidual += residual;

	// handle particles
	if (_comm->getRank() == 0){
		_particles->timestep(dt, _u, _v);
		_particles->writeToFile(_t);
	}

	// print information
	if (printInfo){

		// print running clock, next position
		_clock->nextPos();

		// set curser back, such that the console output it being overwritten
		std::cout << "\r\x1b[A\x1b[A\x1b[A\x1b[A\x1b[A";
		std::cout << "\x1b[A\x1b[A\x1b[A\x1b[A\x1b[A\x1b[A";

		// the actual console output
		std::cout << "============================================================\n";

		std::cout << "Progress: \t\t\t";
		printf("%.2f", _t / _param->Tend() * 100.0);
		std::cout << "\t %" << std::flush;
		std::cout << "\t\t" << _clock->repr(0) << "\n";

		std::cout << "Total simulated time: t = \t";
		printf("%.2f", _t);
		std::cout << "  \t seconds"; // total simulated time
		std::cout << "\t" << _clock->repr(1) << "\n";

		if (_solver_converging){
			std::cout << "(Solver is converging.)        " << std::flush;
		} else {
			std::cout << "(Solver is NOT converging!)    " << std::flush;
		}
		std::cout << "\t\t\t\t" << _clock->repr(2) << "\n";

		std::cout << "Last residual: res = ";
		printf("%10.4f", currentResidual/currentNoTimeSteps);
		std::cout << ",    \tno. iterations: ";
		printf("%7i", int(currentIterations/currentNoTimeSteps/2));
		std::cout << "     \n"; // residual

		std::cout << "Last timestep: dt =  ";
		printf("%10.4f", currentTime/currentNoTimeSteps);
		std::cout << ",    \tno. timesteps:  ";
		printf("%7i", int(currentNoTimeSteps));
		std::cout << "     \n"; // timestep
		currentResidual = 0.0;
		currentNoTimeSteps = 0.0;
		currentIterations = 0.0;
		currentTime = 0.0;

		// magnitudes of the fields
		std::cout << "max(F) = ";
		printf("%7.4f", _F->TotalAbsMax());
		std::cout << ", \tmax(G) = ";
		printf("%7.4f", _G->TotalAbsMax());
		std::cout << ", \tmax(rhs) = ";
		printf("%7.4f", _rhs->TotalAbsMax());
		std::cout << "   \n";
		std::cout << "max(u) = ";
		printf("%7.4f", _u->TotalAbsMax());
		std::cout << ", \tmax(v) = ";
		printf("%7.4f", _v->TotalAbsMax());
		std::cout << ", \tmax(p) = ";
		printf("%7.4f", _p->TotalAbsMax());
		std::cout << "   \n";
		
		//std::cout << "Average value of rhs: " << _rhs->average_value() << "\n";

		std::cout << "diff(F) = ";
		printf("%10.7f", _diff_F);
		std::cout << ", \tdiff(G) = ";
		printf("%10.7f", _diff_G);
		std::cout << std::endl;
		std::cout << "diff(p) = ";
		printf("%10.7f", _diff_p);
		std::cout << ", \tdiff(rhs) = ";
		printf("%10.7f", _diff_rhs);
		std::cout << std::endl;

		std::cout << "============================================================\n";

	} else {
		_F->TotalAbsMax();
		_G->TotalAbsMax();
		_rhs->TotalAbsMax();
		_u->TotalAbsMax();
		_v->TotalAbsMax();
		_p->TotalAbsMax();
	}

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
	//return _rhs; //TODO REMOVE!!!!
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
\return pointer on the Grid containing the stream function evaluations
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
	InteriorIteratorGG it(_geom);
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
	InteriorIteratorGG it(_geom);
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
	InteriorIteratorGG it(_geom);
	it.First();
	while (it.Valid()){
		_rhs->Cell(it) = 1.0/dt * (_F->dx_l(it) + _G->dy_l(it));
		//_rhs->Cell(it) = 1.0/dt * (_F->dx_c(it) + _G->dy_c(it));

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

/**
Update the boundary values in all grids where this is necessary.
*/
void Compute::update_boundary_values()
{
	_geom->UpdateGG_U(_u);
	_geom->UpdateGG_V(_v);

	_geom->UpdateGG_U(_F);
	_geom->UpdateGG_V(_G);
	//_geom->UpdateGG_F(_F,_u);
	//_geom->UpdateGG_G(_G,_v);
}

void Compute::sync_FG()
{
	if (_comm->getSize() > 1) {
		_comm->copyBoundary(_F);
		_comm->copyBoundary(_G);
	}
}

void Compute::sync_uv()
{
	if (_comm->getSize() > 1) {
		_comm->copyBoundary(_u);
		_comm->copyBoundary(_v);
	}
}

void Compute::sync_p()
{
	if (_comm->getSize() > 1) {
		_comm->copyBoundary(_p);
	}
}

void Compute::sync_all()
{
	sync_FG();
	sync_uv();
	sync_p();
}

void Compute::check_for_incontinuities()
{
	real_t diff_p(0.0);
	InteriorIteratorGG it(_geom);
	it.First();
	while (it.Valid()){
		diff_p += std::abs(_p->Cell(it) - _p_old->Cell(it));
		it.Next();
	}
	diff_p /= _geom->Size()[0] * _geom->Size()[1];

	_diff_p = diff_p;
	if (diff_p > 0.2 && _t > 0.1){
//		for (int i=0; i<10; i++) std::cout << "Warning: Incontinuity (P) detected! Diff: " << diff_p << std::endl;
	}

	real_t diff_rhs = 0.0;
	InteriorIteratorGG it2(_geom);
	it2.First();
	while (it2.Valid()){
		diff_rhs += std::abs(_rhs->Cell(it2) - _rhs_old->Cell(it2));
		it2.Next();
	}
	diff_rhs /= _geom->Size()[0] * _geom->Size()[1];

	_diff_rhs = diff_rhs;
	if (diff_rhs > 0.2 && _t > 0.1){
//		for (int i=0; i<10; i++) std::cout << "Warning: Incontinuity (RHS) detected! Diff: " << diff_rhs << std::endl;
	}

	real_t diff_F = 0.0;
	InteriorIteratorGG it3(_geom);
	it3.First();
	while (it3.Valid()){
		diff_F += std::abs(_F->Cell(it3) - _F_old->Cell(it3));
		it3.Next();
	}
	diff_F /= _geom->Size()[0] * _geom->Size()[1];

	_diff_F = diff_F;
	if (diff_F > 0.2 && _t > 0.1){
//		for (int i=0; i<10; i++) std::cout << "Warning: Incontinuity (F) detected! Diff: " << diff_F << std::endl;
	}

	real_t diff_G = 0.0;
	InteriorIteratorGG it4(_geom);
	it4.First();
	while (it4.Valid()){
		diff_G += std::abs(_G->Cell(it4) - _G_old->Cell(it4));
		it4.Next();
	}
	diff_G /= _geom->Size()[0] * _geom->Size()[1];

	_diff_G = diff_G;
	if (diff_G > 0.2 && _t > 0.1){
//		for (int i=0; i<10; i++) std::cout << "Warning: Incontinuity (G) detected! Diff: " << diff_G << std::endl;
	}

}
