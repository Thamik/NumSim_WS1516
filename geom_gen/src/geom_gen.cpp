#include "geom_gen.hpp"

#include <iostream>
#include <fstream>

/**
The default constructor of the GeometryGenerator class.
*/
GeometryGenerator::GeometryGenerator()
: _bSizeX(128), _bSizeY(128), _bLengthX(1.0), _bLengthY(1.0), _filename(nullptr), _flags(nullptr), _bvu(nullptr), _bvv(nullptr), _bvp(nullptr)
{
	setSize(_bSizeX, _bSizeY);
	setLength(_bLengthX, _bLengthY);
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

void GeometryGenerator::initZero()
{
	for (int i=0; i<_bSizeX*_bSizeY; i++){
		_flags[i] = 0;
		_bvu[i] = 0;
		_bvv[i] = 0;
		_bvp[i] = 0;
	}
}

void GeometryGenerator::autoBalance()
{
	// TODO
}

void GeometryGenerator::drivenCavity()
{
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
	// assure that all values are initialized with zero
	initZero();

	setLength(xlength, ylength);
	
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
	int stepSizeX = int(_bSizeX * stepLength/_bLengthX);
	int stepSizeY = int(_bSizeY * stepLength/_bLengthY);

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
}

void karmanVortexStreet(double alpha, double width, double xlength, double ylength, double pressureLeft, double pressureRight)
{
	// assume that ylength is greater or equal to xlength
	if (ylength<xlength){
		std::cout << "Warning: Karman Vortex Street geometry could not be constructed. ylength < xlength!" << std::endl;
	}

	pipeFlow(xlength, ylength, pressureLeft, pressureRight);

	double length/2.0;
	
}
