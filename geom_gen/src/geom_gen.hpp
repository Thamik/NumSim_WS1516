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

	void drivenCavity();

	void writeToFile(const char* filename=nullptr);
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
