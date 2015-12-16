#include <string>

/**
This class implements an complex geometry generator.
*/
class GeometryGenerator {
public:
	GeometryGenerator();
	~GeometryGenerator();

	void setTotalSize(int n);
	void setLength(double x, double y);

	void writeToFile(const char* filename=nullptr, const char* filenameParam=nullptr);

	void print() const;

	void fixSingleCells();

	/// Constructs geometry data for the driven cavity problem
	void drivenCavity();
	/// Constructs geometry data for a simple pipe flow
	void pipeFlow(double xlength=6.0, double ylength=1.0, double pressureLeft=1.1, double pressureRight=1.0);
	/// Constructs geometry data for the flow over a staircase in a pipe
	void flowOverAStep(double xlength=6.0, double ylength=1.0, double pressureLeft=1.1, double pressureRight=1.0);

	void testCase1();

	void testCase2();

	void karmanVortexStreet(double alpha, double width=0.2, double xlength=6.0, double ylength=1.0, double pressureLeft=1.1, double pressureRight=1.0);

private:
	int _bSizeX, _bSizeY;
	int _totalCells;
	double _bLengthX, _bLengthY;

	std::string* _filename;
	std::string* _filenameParam;

	char* _flags;
	double* _bvu;
	double* _bvv;
	double* _bvp;

	// local Parameter variables
	double _re; // Reynolds number
  	double _omega; // relaxation factor
  	double _alpha; // upwind differencing factor
  	double _dt; // time step size
  	double _tend; // end time
  	double _eps; // tolerance for pressure iteration
  	double _tau; // safety factor time step size
  	int _itermax; // maximum number of iterations
	int _itermin; // minimum number of iterations

	void setSize(int x, int y);

	void autoBalance();

	void initZero();

	bool isObstacle(int ival) const;
	bool isObstacle(int x, int y) const;

	void setNoSlip(int x, int y);
	void set(int x, int y, bool obstacle, bool condu, bool condv, bool condp, double valu, double valv, double valp);
	void set(int x, int y, char flag, double valu, double valv, double valp);

	// help function for the computation of the Kalman vortex street
	bool isInsideObstacle(double alpha, double width, double length, double x0, double y0, double x, double y) const;
};
