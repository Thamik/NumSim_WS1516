#include "typedef.hpp"

/**
This class implements a doubly-linked list, in which each node stores the data concerning a single particle.
*/
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

	bool isInside(multi_real_t pos) const;

private:
	multi_real_t _pos;

	real_t _radius;

	ParticleList* _pred;
	ParticleList* _succ;

	const Geometry* _geom;

};

//--------------------------------------------------------------------

/**
This class stores the data concerning all particles that are involved in the current simulation.
*/
class Particles {
public:
	Particles(const Geometry* geom);
	~Particles();

	void updateAllPositions(real_t dt, const Grid* u, const Grid* v);

	void spawnParticle(multi_real_t pos);

	void streaklinePolicy();
	void particleTracingPolicy();

	void setDefaultFormat();
	void setMatlabFormat();
	void setMatlabOneFileFormat();
	void setPythonOneFileFormat();

	void init();

	void timestep(real_t dt, const Grid* u, const Grid* v);

	void newParticles();

	bool isInsideParticle(multi_real_t pos) const;

	void writeToFile(real_t total_time = -1.0, const char* filename = "") const;
	void finalizeFile() const;

	bool isEmpty() const;

	int numberParticles() const;

private:
	ParticleList* _pl;

	const Geometry* _geom;

	enum {
		undefined = 0,
		streakline = 1,
		particleTracing = 2
	};
	char _policy;

	enum {
		defaultFormat = 0,
		matlabFormat = 1,
		matlabOneFileFormat = 2,
		pythonOneFileFormat = 3
	};
	char _file_format;

	int _no_timestep;

	void init_particleTracing();
	void init_streakline();

	int _no_streaklines;
	int _max_particles_per_streakline;
	multi_real_t* _streakline_positions;

	int _max_particles;

};

