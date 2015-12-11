#include "particles.hpp"
#include "grid.hpp"

/**
Constructor
*/
ParticleList::ParticleList(multi_real_t pos, const Geometry* geom)
: _pos(pos), _pred(nullptr), _succ(nullptr), _geom(geom)
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

//--------------------------------------------------------------------

Particles::Particles(const Geometry* geom)
: _pl(nullptr), _geom(geom), _policy(0)
{
}

Particles::~Particles()
{
	if (_pl != nullptr){
		_pl->deleteAllAfter();
		_pl = nullptr;
	}
}

void Particles::updateAllPositions(real_t dt, const Grid* u, const Grid* v)
{
	ParticleList* particle = _pl;
	while (particle != nullptr){
		particle->updatePos(dt,u,v);
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

void Particles::timestep()
{
	
}

