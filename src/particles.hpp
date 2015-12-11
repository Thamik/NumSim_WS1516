#include "typedef.hpp"

class ParticleList {
public:
	ParticleList(multi_real_t pos, const Geometry* geom);
	~ParticleList();

	multi_real_t& getPos();
	const multi_real_t& getPos() const;

	bool isSingleton() const;
	bool isEnd() const;
	bool isBegin() const;

	void append(multi_real_t pos);

	void deleteSingleParticle();
	void deleteAllAfter();

	void setPred(ParticleList* pred);
	void setSucc(ParticleList* succ);

	ParticleList* getPred();
	ParticleList* getSucc();

	void updatePos(real_t dt, const Grid* u, const Grid* v);

private:
	multi_real_t _pos;

	ParticleList* _pred;
	ParticleList* _succ;

	const Geometry* _geom;

};

//--------------------------------------------------------------------

class Particles {
public:
	Particles(const Geometry* geom);
	~Particles();

	void updateAllPositions(real_t dt, const Grid* u, const Grid* v);

	void spawnParticle(multi_real_t pos);

	void streaklinePolicy();
	void particleTracingPolicy();

	void timestep();

private:
	ParticleList* _pl;

	const Geometry* _geom;

	enum {
		undefined = 0,
		streakline = 1,
		particleTracing = 2
	};
	char _policy;

};

