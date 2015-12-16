#include <iostream>
#include <fstream>
#include <math.h>
#include <string>

#include "particles.hpp"
#include "grid.hpp"
#include "geometry.hpp"

/**
Constructor
*/
ParticleList::ParticleList(multi_real_t pos, const Geometry* geom)
: _pos(pos), _radius(0.02), _pred(nullptr), _succ(nullptr), _geom(geom)
{
}

/**
Destructor
*/
ParticleList::~ParticleList()
{
}

bool ParticleList::isSingleton() const
{
	return isEnd() && isBegin();
}

bool ParticleList::isEnd() const
{
	return _succ == nullptr;
}

bool ParticleList::isBegin() const
{
	return _pred == nullptr;
}

void ParticleList::append(multi_real_t pos)
{
	if (isEnd()){
		_succ = new ParticleList(pos,_geom);
		_succ->setPred(this);
	} else {
		_succ->append(pos);
	}
}

void ParticleList::setPred(ParticleList* pred)
{
	_pred = pred;
}

void ParticleList::setSucc(ParticleList* succ)
{
	_succ = succ;
}

void ParticleList::deleteSingleParticle()
{
	ParticleList* temp_pred = _pred;
	ParticleList* temp_succ = _succ;

	if (!isBegin()){
		temp_pred->setSucc(temp_succ);
	}
	if (!isEnd()){
		temp_succ->setPred(temp_pred);
	}

	delete this; // delete self
}

void ParticleList::deleteAllAfter()
{
	if (!isBegin()){
		_pred->setSucc(nullptr);
	}

	if (!isEnd()){
		_succ->deleteAllAfter();
	}

	delete this;
}

multi_real_t& ParticleList::getPos()
{
	return _pos;
}

const multi_real_t& ParticleList::getPos() const
{
	return _pos;
}

ParticleList* ParticleList::getPred()
{
	return _pred;
}

ParticleList* ParticleList::getSucc()
{
	return _succ;
}

void ParticleList::updatePos(real_t dt, const Grid* u, const Grid* v)
{
	multi_real_t temp_pos = _pos;
	_pos[0] += dt * u->Interpolate(temp_pos);
	_pos[1] += dt * v->Interpolate(temp_pos);
}

bool ParticleList::isInside(multi_real_t pos) const
{
	return sqrt(pow(pos[0]-_pos[0], 2.0) - pow(pos[1]-_pos[1], 2.0)) <= _radius;
}

//--------------------------------------------------------------------

Particles::Particles(const Geometry* geom)
: _pl(nullptr), _geom(geom), _policy(0), _file_format(defaultFormat), _no_timestep(0), _no_streaklines(0), _max_particles_per_streakline(0), _streakline_positions(nullptr), _max_particles(0)
{
}

Particles::~Particles()
{
	if (_pl != nullptr){
		_pl->deleteAllAfter();
		_pl = nullptr;
	}

	if (_streakline_positions != nullptr) delete[] _streakline_positions;
}

void Particles::updateAllPositions(real_t dt, const Grid* u, const Grid* v)
{
	ParticleList* particle = _pl;
	while (particle != nullptr){
		if (_geom->isInsideThisSubdomain(particle->getPos())){
			particle->updatePos(dt,u,v);
		}
		particle = particle->getSucc();
	}
}

void Particles::spawnParticle(multi_real_t pos)
{
	if (_pl == nullptr){
		_pl = new ParticleList(pos,_geom);
	} else {
		_pl->append(pos);
	}
}

void Particles::streaklinePolicy()
{
	_policy = streakline;
}

void Particles::particleTracingPolicy()
{
	_policy = particleTracing;
}

void Particles::init()
{
	switch (_policy){
		case streakline:
			init_streakline();
			break;
		case particleTracing:
			init_particleTracing();
			break;
		case undefined:
		default:
			std::cout << "Warning! Particles: Unknown or undefined policy! Abort init." << std::endl;
			return;
	}

	if (_file_format == matlabOneFileFormat){
		remove("VTK/timestep_all_particles.m");
	} else if (_file_format == pythonOneFileFormat){
		remove("VTK/timestep_all_particles.py");

		std::string file("");
		file += "VTK/timestep_all_particles.py";
		std::ofstream outfile;
		outfile.open(file.c_str());

		outfile << "import numpy as np\n";
		outfile << "pos = np.array([";

		outfile.close();
	}
}

void Particles::init_streakline()
{
	if (_streakline_positions != nullptr) delete[] _streakline_positions;

	_no_streaklines = 4;
	_max_particles_per_streakline = 500;
	_streakline_positions = new multi_real_t[_no_streaklines];

	// in the middle
//	_streakline_positions[0] = multi_real_t(_geom->TotalLength()[0]/2.0, _geom->TotalLength()[1]/2.0);

	// in the left quarter, which is very nice for pipe flow (with or without step) and karman vortex street
	//_streakline_positions[0] = multi_real_t(_geom->TotalLength()[0]/4.0, _geom->TotalLength()[1]/2.0);
	_streakline_positions[0] = multi_real_t(0.2, 0.8);
	_streakline_positions[1] = multi_real_t(0.2, 0.6);
	_streakline_positions[2] = multi_real_t(0.2, 0.4);
	_streakline_positions[3] = multi_real_t(0.2, 0.2);
	

	_max_particles = _no_streaklines * _max_particles_per_streakline;
}

void Particles::init_particleTracing()
{
	// spawn all particles
	
	// for instance, spawn the particles on a uniform grid all over the geometry
	int totalNumberParticles(400);
	int nx, ny;
	// split the number of particles in each dimension like the mesh was built
	nx = int(round(real_t(_geom->TotalSize()[0]) * sqrt(real_t(totalNumberParticles) / real_t(_geom->TotalSize()[0] * _geom->TotalSize()[1]))));
	ny = int(round(real_t(_geom->TotalSize()[1]) * sqrt(real_t(totalNumberParticles) / real_t(_geom->TotalSize()[0] * _geom->TotalSize()[1]))));

	for (int i=1; i<=nx; i++){
		for (int j=1; j<=ny; j++){
			real_t x = _geom->TotalLength()[0]/(nx+1.0)*real_t(i);
			real_t y = _geom->TotalLength()[1]/(ny+1.0)*real_t(j);
			spawnParticle(multi_real_t(x,y));
		}
	}
	std::cout << "Particle Tracing: All " << numberParticles() << " particles spawned\n";

	_max_particles = totalNumberParticles;
}

void Particles::timestep(real_t dt, const Grid* u, const Grid* v)
{
	_no_timestep++;

	updateAllPositions(dt, u, v);
	newParticles();
}

void Particles::newParticles()
{
	switch (_policy){
		case streakline:
			// spawn particles at the specified positions
			for (int i=0; i<_no_streaklines; i++){
				if (numberParticles() < _max_particles){
					spawnParticle(_streakline_positions[i]);
				} else {
					break;
				}
			}
			break;
		case particleTracing:
			// nothing to do here, all particles spawned in the beginning
			break;
		case undefined:
		default:
			std::cout << "Warning! Particles: Unknown or undefined policy!" << std::endl;
	}
}

bool Particles::isInsideParticle(multi_real_t pos) const
{
	ParticleList* particle = _pl;
	while (particle != nullptr){
		if (particle->isInside(pos)) return true;
		particle = particle->getSucc();
	}
	return false;
}

void Particles::writeToFile(real_t total_time, const char* filename) const
{
	std::string file("");
	if (_file_format == defaultFormat || _file_format == matlabFormat){
		if (!strcmp(filename, "")){
			file += "VTK/timestep_";
			file += std::to_string(_no_timestep);
			file += ".particles";
		} else {
			file += filename;
		}
	} else if (_file_format == matlabOneFileFormat){
		file += "VTK/timestep_all_particles.m";
	} else if (_file_format == pythonOneFileFormat){
		file += "VTK/timestep_all_particles.py";
	} else {
		std::cout << "Warning! Particles: Unknown file format. Abort writing." << std::endl;
		return;
	}

	//std::cout << "Writing particle data to file: " << file << std::endl;
	//std::cout << "In total " << numberParticles() << " particles." << std::endl;

	std::ofstream outfile;
	if (_file_format == defaultFormat || _file_format == matlabFormat){
		outfile.open(file.c_str());
	} else if (_file_format == matlabOneFileFormat || _file_format == pythonOneFileFormat){
		outfile.open(file.c_str(), std::ofstream::app);
	}

	// check on filestream
	if (!outfile.is_open()){
		// something went wrong, file is not open
		std::cout << "Warning! Particles: Could not open file! No output file was written!\n" << std::flush;
		return;
	}

	// write title
	if (_file_format == defaultFormat){
		outfile << "Particle Positions\n";
	}
	// write total time
	if (_file_format == defaultFormat){
		outfile << "Total time: t = ";
		if (total_time < 0.0){
			outfile << "-";
		} else {
			outfile << total_time;
		}
		outfile << "\n";
		// one empty line
		outfile << "\n";
	}
	// write particle positions
	if (_file_format == matlabFormat || _file_format == matlabOneFileFormat){
		outfile << "pos(:,:," << _no_timestep << ") = [ ";
	} else if (_file_format == pythonOneFileFormat){
		outfile << "[";
	}

	ParticleList* particle = _pl;
	while (particle != nullptr){
		if (_file_format == defaultFormat){
			outfile << particle->getPos()[0] << " " << particle->getPos()[1] << "\n";
		} else if (_file_format == matlabFormat || _file_format == matlabOneFileFormat){
			outfile << particle->getPos()[0] << " , " << particle->getPos()[1];
			if (!particle->isEnd()){
				outfile << " ;\n";
			}
		} else if (_file_format == pythonOneFileFormat){
			outfile << "[ " << particle->getPos()[0] << ", " << particle->getPos()[1] << " ],\n";
		}

		particle = particle->getSucc();
	}
	if (_file_format == matlabFormat || _file_format == matlabOneFileFormat){
		outfile << " ];";
	}
	if (_file_format == matlabOneFileFormat){
		outfile << "\n\n";
	}
	if (_file_format == pythonOneFileFormat){
		if (numberParticles() < _max_particles){
			// write "empty" particle data, such that the (sub-)array has the right shape
			int diff_particles = _max_particles - numberParticles();
			for (int i=0; i<diff_particles; i++){
				outfile << "[ float(\"NaN\"), float(\"NaN\") ],\n";
			}
		} else if (numberParticles() > _max_particles){
			// something went wrong, there are too much particles!
			std::cout << "Warning! Particles: too much particles." << std::endl;
		}
		outfile << "],\n";
	}

	outfile.close();
}

void Particles::finalizeFile() const
{
	if (_file_format == pythonOneFileFormat){
		std::string file("");
		file += "VTK/timestep_all_particles.py";
		std::ofstream outfile;
		outfile.open(file.c_str(), std::ofstream::app);

		outfile << "])";

		outfile.close();
	}
}

bool Particles::isEmpty() const
{
	return _pl == nullptr;
}

int Particles::numberParticles() const
{
	int res(0);
	ParticleList* particle = _pl;
	while (particle != nullptr){
		res++;
		particle = particle->getSucc();
	}
	return res;
}

void Particles::setDefaultFormat()
{
	_file_format = defaultFormat;
}

void Particles::setMatlabFormat()
{
	_file_format = matlabFormat;
}

void Particles::setMatlabOneFileFormat()
{
	_file_format = matlabOneFileFormat;
}

void Particles::setPythonOneFileFormat()
{
	_file_format = pythonOneFileFormat;
}
