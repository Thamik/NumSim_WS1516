#include <string>

/**
This class implements an complex geometry generator.
*/
class GeometryGenerator {
public:
	GeometryGenerator();
	~GeometryGenerator();

	void setSize(int x, int y);
	void setLength(double x, double y);

	void writeToFile(const char* filename=nullptr);

	void autoBalance();

	/// Constructs geometry data for the driven cavity problem
	void drivenCavity();
	/// Constructs geometry data for a simple pipe flow
	void pipeFlow(double xlength=6.0, double ylength=1.0, double pressureLeft=2, double pressureRight=1);
	/// Constructs geometry data for the flow over a staircase in a pipe
	void flowOverAStep(double xlength=6.0, double ylength=1.0, double pressureLeft=2, double pressureRight=1);

	void karmanVortexStreet(double alpha, double width=0.2, double xlength=6.0, double ylength=1.0, double pressureLeft=2, double pressureRight=1);

private:
	int _bSizeX, _bSizeY;
	double _bLengthX, _bLengthY;

	std::string* _filename;

	char* _flags;
	double* _bvu;
	double* _bvv;
	double* _bvp;

	void initZero();
};
