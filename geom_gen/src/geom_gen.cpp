#include "geom_gen.hpp"

#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm> // std::min, std::max

/**
The default constructor of the GeometryGenerator class.
*/
GeometryGenerator::GeometryGenerator()
: _bSizeX(1), _bSizeY(1), _totalCells(10000), _bLengthX(1.0), _bLengthY(1.0), _filename(nullptr), _flags(nullptr), _bvu(nullptr), _bvv(nullptr), _bvp(nullptr)
{
	setLength(_bLengthX, _bLengthY);
	setTotalSize(_totalCells);
	drivenCavity();
}

/**
The destructor of the GeometryGenerator class.
*/
GeometryGenerator::~GeometryGenerator()
{
	if (_filename != nullptr) delete _filename;
	if (_flags != nullptr) delete _flags;
	if (_bvu != nullptr) delete _bvu;
	if (_bvv != nullptr) delete _bvv;
	if (_bvp != nullptr) delete _bvp;
}

/**
Sets the total sizes of the goemetry.
*/
void GeometryGenerator::setSize(int x, int y)
{
	_bSizeX = x+2;
	_bSizeY = y+2;
	_totalCells = _bSizeX *_bSizeY;
	// delete memory if necessary
	if (_flags != nullptr) delete _flags;
	if (_bvu != nullptr) delete _bvu;
	if (_bvv != nullptr) delete _bvv;
	if (_bvp != nullptr) delete _bvp;
	// allocate memory for all the data fields
	_flags = new char[_bSizeX * _bSizeY];
	_bvu = new double[_bSizeX * _bSizeY];
	_bvv = new double[_bSizeX * _bSizeY];
	_bvp = new double[_bSizeX * _bSizeY];
	// initialize fields with zero
	initZero();
}

void GeometryGenerator::setTotalSize(int n)
{
	_totalCells = n;
	autoBalance();
}

void GeometryGenerator::setLength(double x, double y)
{
	_bLengthX = x;
	_bLengthY = y;
}

/**
Writes the complex geometry data to a specified file.
*/
void GeometryGenerator::writeToFile(const char* filename)
{
	if (filename == nullptr && _filename == nullptr){
		// handle default filename parameter, set this value to default
//		_filename = new std::string("../data/complex_default.geom");
		_filename = new std::string("geom_files/complex_default.geom");
		std::cout << "Writing geometry data to (default) file: " << _filename->c_str() << std::endl;
	} else if (filename == nullptr) {
		// nothing given, but already a name in _filename (member)
		std::cout << "Writing geometry data to file: " << _filename->c_str() << std::endl;
		// just use the path in _filename, do nothing
	} else {
		if (_filename != nullptr) delete _filename;
		_filename = new std::string(filename);
		std::cout << "Writing geometry data to file: " << _filename->c_str() << std::endl;
	}

	// write to file
	std::ofstream outfile;
	outfile.open(_filename->c_str());
	if (!outfile.is_open()){
		// something went wrong, file is not open
		std::cout << "Warning: Geometry generator could not open file! No output file was written!\n" << std::flush;
	} else {
		// file is open, go on writing
		outfile << _bSizeX-2 << "\n"; // write without ghost cells
		outfile << _bSizeY-2 << "\n";
		outfile << _bLengthX << "\n";
		outfile << _bLengthY << "\n";
		
		for (int i=0; i<_bSizeX*_bSizeY; i++){
			outfile << int(_flags[i]) << " ";
		}
		outfile << "\n";
		for (int i=0; i<_bSizeX*_bSizeY; i++){
			outfile << _bvu[i] << " ";
		}
		outfile << "\n";
		for (int i=0; i<_bSizeX*_bSizeY; i++){
			outfile << _bvv[i] << " ";
		}
		outfile << "\n";
		for (int i=0; i<_bSizeX*_bSizeY; i++){
			outfile << _bvp[i] << " ";
		}
//		outfile << "\n";
	}
	outfile.close();
}

void GeometryGenerator::print() const
{
	std::cout << "==================================================\n";
	std::cout << "GeomGen: Geometry configuration:\n";
	std::cout << "Total size, without ghost cells: (" << _bSizeX-2 << ", " << _bSizeY-2 << ")\n";
	std::cout << "Total length: (" << _bLengthX << ", " << _bLengthY << ")\n";
	std::cout << "==================================================\n" << std::flush;
}

void GeometryGenerator::initZero()
{
	for (int i=0; i<_bSizeX*_bSizeY; i++){
		_flags[i] = 0;
		_bvu[i] = 0;
		_bvv[i] = 0;
		_bvp[i] = 0;
	}
}

bool GeometryGenerator::isObstacle(int ival) const
{
	return _flags[ival] & 1;
}

bool GeometryGenerator::isObstacle(int x, int y) const
{
	return isObstacle(y * (_bSizeX) + x);
}

void GeometryGenerator::setNoSlip(int x, int y)
{
	set(x,y,1,0,0,1,0.0,0.0,0.0);
}

void GeometryGenerator::set(int x, int y, bool obstacle, bool condu, bool condv, bool condp, double valu, double valv, double valp)
{
	char flag = char(obstacle) | (char(condu) << 1) | (char(condv) << 2) | (char(condp) << 3);

	set(x,y,flag,valu,valv,valp);
}

void GeometryGenerator::set(int x, int y, char flag, double valu, double valv, double valp)
{
	int ival = y * (_bSizeX) + x;
	_flags[ival] = flag;
	_bvu[ival] = valu;
	_bvv[ival] = valv;
	_bvp[ival] = valp;
}

void GeometryGenerator::autoBalance()
{
	// try to choose the size in x- and y-direction such that the problem
	// is more balanced and the total number of cells is being conserved
	if (_totalCells <= 0) _totalCells = _bSizeX * _bSizeY;

	_bSizeX = int(round( sqrt(_bLengthX/_bLengthY * _totalCells) ));
	_bSizeY = int(round( sqrt(_bLengthY/_bLengthX * _totalCells) ));

	setSize(_bSizeX, _bSizeY);
}

void GeometryGenerator::fixSingleCells()
{
	// iterate over all interior cells and check for problems
	int ival(0);
	for (int i=1; i<_bSizeX-1; i++){
		for (int j=1; j<_bSizeY-1; j++){
			ival = j * (_bSizeX) + i;
			
			if (isObstacle(i,j) && ((!isObstacle(i-1,j) && !isObstacle(i+1,j)) || (!isObstacle(i,j-1) && !isObstacle(i,j+1)))){
				// critical cell

				// fix: add 3 points in L-form around this cell
				char flag = _flags[ival];
				double valu = _bvu[ival];
				double valv = _bvv[ival];
				double valp = _bvp[ival];
				// right cell
				set(i+1,j,flag,valu,valv,valp);
				// lower cell
				set(i,j-1,flag,valu,valv,valp);
				// lower right cell
				set(i+1,j-1,flag,valu,valv,valp);
			}
		}
	}
}

void GeometryGenerator::drivenCavity()
{
	// auto balance the number of cells in each dimension
	autoBalance();

	// assure that all values are initialized with zero
	initZero();

	int ival(0);
	for (int i=0; i<_bSizeX; i++){
		for (int j=0; j<_bSizeY; j++){
			ival = j * (_bSizeX) + i;
			if (j==_bSizeY-1) {
				// upper boundary
				_flags[ival] = 1 | 1<<3; // neumann condition for p, dirichlet for u and v
				_bvu[ival] = 1.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = 0.0;
			} else if (j==0 || i==0 || i==_bSizeX-1){
				// lower, left or right boundary
				_flags[ival] = 1 | 1<<3; // neumann condition for p, dirichlet for u and v
				_bvu[ival] = 0.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = 0.0;
			}
		}
	}
}

void GeometryGenerator::pipeFlow(double xlength, double ylength, double pressureLeft, double pressureRight)
{
	setLength(xlength, ylength);

	// auto balance the number of cells in each dimension
	autoBalance();

	// assure that all values are initialized with zero
	initZero();
	
	int ival(0);
	for (int i=0; i<_bSizeX; i++){
		for (int j=0; j<_bSizeY; j++){
			ival = j * (_bSizeX) + i;
			if (j==_bSizeY-1) {
				// upper boundary
				_flags[ival] = 1 | 1<<3; // neumann condition for p, dirichlet for u and v
				_bvu[ival] = 0.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = 0.0;
			} else if (j==0){
				// lower boundary
				_flags[ival] = 1 | 1<<3; // neumann condition for p, dirichlet for u and v
				_bvu[ival] = 0.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = 0.0;
			} else if (i==0){
				// left boundary
				_flags[ival] = 1 | 1<<1 | 1<<2; // dirichlet condition for p, neumann for u and v
				_bvu[ival] = 0.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = pressureLeft;
			} else if (i==_bSizeX-1){
				// right boundary
				_flags[ival] = 1 | 1<<1 | 1<<2; // dirichlet condition for p, neumann for u and v
				_bvu[ival] = 0.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = pressureRight;
			}
		}
	}
}

void GeometryGenerator::flowOverAStep(double xlength, double ylength, double pressureLeft, double pressureRight)
{
	pipeFlow(xlength, ylength, pressureLeft, pressureRight);
	
	double stepLength = ylength/2.0;
	int stepSizeX = int(round(_bSizeX * stepLength/_bLengthX));
	int stepSizeY = int(round(_bSizeY * stepLength/_bLengthY));

	// set all cells that describe the step
	int ival(0);
	for (int i=0; i<stepSizeX; i++){
		for (int j=0; j<stepSizeY; j++){
			ival = j * (_bSizeX) + i;

			// set NO-SLIP conditions
			_flags[ival] = 1 | 1<<3; // neumann condition for p, dirichlet for u and v
			_bvu[ival] = 0.0;
			_bvv[ival] = 0.0;
			_bvp[ival] = 0.0;
		}
	}

	std::cout << pressureLeft << ", " << pressureRight << "\n" << std::flush;
}

void GeometryGenerator::karmanVortexStreet(double alpha, double width, double xlength, double ylength, double pressureLeft, double pressureRight)
{
	// assume that ylength is greater or equal to xlength
	if (xlength<ylength){
		std::cout << "Warning: Karman Vortex Street geometry could not be constructed. xlength < ylength!" << std::endl;
		return;
	}

	pipeFlow(xlength, ylength, pressureLeft, pressureRight);

	double length = ylength/2.0;

	// set all cells that describe the obstacle that causes the vortex street
	int cellsX = int(round(2.0*length/xlength * _bSizeX));
	int cellsY = int(round(2.0*length/ylength * _bSizeY));
	cellsX = std::min(cellsX, _bSizeX);
	cellsY = std::min(cellsY, _bSizeY);

	int ival(0);
	double xPos(0.0);
	double yPos(0.0);
	for (int i=0; i<cellsX; i++){
		for (int j=0; j<cellsY; j++){
			ival = j * (_bSizeX) + i;
			xPos = double(i)/double(_bSizeX-1) * xlength;
			yPos = double(j)/double(_bSizeY-1) * ylength;

			if (isInsideObstacle(alpha, width, length, length, length, xPos, yPos)){
				// set NO-SLIP conditions
				_flags[ival] = 1 | 1<<3; // neumann condition for p, dirichlet for u and v
				_bvu[ival] = 0.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = 0.0;
			}

		}
	}

	fixSingleCells();
	
}

bool GeometryGenerator::isInsideObstacle(double alpha, double width, double length, double x0, double y0, double x, double y) const
{
	double a1 = sin(alpha);
	double a2 = -cos(alpha);
	double b1 = cos(alpha);
	double b2 = sin(alpha);
	double c = x0 - x;
	double d = y0 - y;

	// solve linear equation system (two-dimensional)
	double s = 1/(a1*b2-b1*a2) * (-a2*c + a1*d);
	//double t = 1/(a1*b2-b1*a2) * (-b2*c + b1*d);

	double proj1 = s*b1 + x;
	double proj2 = s*b2 + y;

	double distance_length = sqrt( pow(proj1-x,2.0) + pow(proj2-y,2.0) );
	double distance_width = sqrt( pow(proj1-x0,2.0) + pow(proj2-y0,2.0) );

	return (distance_length <= length/2.0) && (distance_width <= width/2.0);
}



void GeometryGenerator::testCase1()
{
	double xlength = 6.0;
	double ylength = 1.0;
	double alpha = 0.0;
	double pressureRight = 0.0;
	double pressureLeft = 1.0;
	double width = 0.5;
	
	// assume that ylength is greater or equal to xlength
	if (xlength<ylength){
		std::cout << "Warning: Karman Vortex Street geometry could not be constructed. xlength < ylength!" << std::endl;
		return;
	}

	pipeFlow(xlength, ylength, pressureLeft, pressureRight);

	double length = ylength/2.0;

	// set all cells that describe the obstacle that causes the vortex street
	int cellsX = int(round(2.0*length/xlength * _bSizeX));
	int cellsY = int(round(2.0*length/ylength * _bSizeY));
	cellsX = std::min(cellsX, _bSizeX);
	cellsY = std::min(cellsY, _bSizeY);

	int ival(0);
	double xPos(0.0);
	double yPos(0.0);
	for (int i=0; i<_bSizeX; i++){
		for (int j=0; j<cellsY; j++){
			ival = j * (_bSizeX) + i;
			xPos = double(i)/double(_bSizeX-1) * xlength;
			yPos = double(j)/double(_bSizeY-1) * ylength;

			if (isInsideObstacle(alpha, width, length, length+2.0, length, xPos, yPos)){
				// set NO-SLIP conditions
				_flags[ival] = 1 | 1<<3; // neumann condition for p, dirichlet for u and v
				_bvu[ival] = 0.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = 0.0;
			}

		}
	}

	fixSingleCells();
	
}
