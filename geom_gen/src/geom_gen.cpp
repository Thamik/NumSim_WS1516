#include "geom_gen.hpp"

#include <iostream>
#include <fstream>

#define _USE_MATH_DEFINES // to get pi via M_PI
#include <math.h>

#include <algorithm> // std::min, std::max

/**
The default constructor of the GeometryGenerator class.
*/
GeometryGenerator::GeometryGenerator()
: _bSizeX(1), _bSizeY(1), _totalCells(10000), _bLengthX(1.0), _bLengthY(1.0), _forcePowerOfTwo(false), _filenameGeom(), _filenameParam(), _flags(nullptr), _bvu(nullptr), _bvv(nullptr), _bvp(nullptr), _bvt(nullptr), _re(1000.0), _omega(1.7), _alpha(0.9), _eps(0.001), _tau(0.9), _itermax(1000), _itermin(0), _dt(0.5), _tend(30.0), _gravityX(0.0), _gravityY(0.0), _pr(3.0), _beta(0.5), _gamma(0.9), _q(0.0)
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
	if (_flags != nullptr) delete[] _flags;
	if (_bvu != nullptr) delete[] _bvu;
	if (_bvv != nullptr) delete[] _bvv;
	if (_bvp != nullptr) delete[] _bvp;
	if (_bvt != nullptr) delete[] _bvt;
}

/**
Sets the total sizes of the goemetry.
*/
void GeometryGenerator::setSize(int x, int y)
{
	// if _forcePowerOfTwo is true, make sure to have powers of two as the number of cells in each dimension
	if (_forcePowerOfTwo){
		_bSizeX = int(round(pow(2,round(log2(double(x))))));
		_bSizeY = int(round(pow(2,round(log2(double(y))))));
	} else {
		_bSizeX = x;
		_bSizeY = y;
	}
	_bSizeX += 2; // add ghost cells
	_bSizeY += 2;
	_totalCells = _bSizeX *_bSizeY;
	// delete memory if necessary
	if (_flags != nullptr) delete[] _flags;
	if (_bvu != nullptr) delete[] _bvu;
	if (_bvv != nullptr) delete[] _bvv;
	if (_bvp != nullptr) delete[] _bvp;
	if (_bvt != nullptr) delete[] _bvt;
	// allocate memory for all the data fields
	_flags = new char[_bSizeX * _bSizeY];
	_bvu = new double[_bSizeX * _bSizeY];
	_bvv = new double[_bSizeX * _bSizeY];
	_bvp = new double[_bSizeX * _bSizeY];
	_bvt = new double[_bSizeX * _bSizeY];
	// initialize fields with zero
	initZero();
}

void GeometryGenerator::setTotalSize(int n)
{
	_totalCells = n;
	//autoBalance();
}

void GeometryGenerator::forcePowerOfTwo()
{
	_forcePowerOfTwo = true;
}

void GeometryGenerator::doNotForcePowerOfTwo()
{
	_forcePowerOfTwo = false;
}

void GeometryGenerator::setLength(double x, double y)
{
	_bLengthX = x;
	_bLengthY = y;
}

/**
Writes the complex geometry data to a specified file.
*/
void GeometryGenerator::writeToFile() const
{
	// write geometry file
	if (_filenameGeom.compare("") == 0){
		// no filename specified, dont write file
		std::cout << "Filename not specified, no geometry file was written.\n";
	} else {
		// filename is specified, write to file
		std::cout << "Writing geometry data to file: " << _filenameGeom.c_str() << std::endl;

		// write to file
		std::ofstream outfile;
		outfile.open(_filenameGeom.c_str());
		if (!outfile.is_open()){
			// something went wrong, file is not open
			std::cout << "Warning: Geometry generator could not open file! No geometry file was written!\n" << std::flush;
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
			outfile << "\n";
			for (int i=0; i<_bSizeX*_bSizeY; i++){
				outfile << _bvt[i] << " ";
			}
		}
		outfile.close();
	}

	// write parameter file
	if (_filenameParam.compare("") == 0){
		// no filename specified, dont write file
		std::cout << "Filename not specified, no parameter file was written.\n";
	} else {
		// filename is specified, write to file
		std::cout << "Writing parameter data to file: " << _filenameParam.c_str() << std::endl;

		// write to file
		std::ofstream outfileParam;
		outfileParam.open(_filenameParam.c_str());
		if (!outfileParam.is_open()){
			// something went wrong, file is not open
			std::cout << "Warning: Geometry generator could not open file! No parameter file was written!\n" << std::flush;
		} else {
			outfileParam << _re << "\n";
			outfileParam << _omega << "\n";
			outfileParam << _alpha << "\n";
			outfileParam << _eps << "\n";
			outfileParam << _tau << "\n";
			outfileParam << _itermax << "\n";
			outfileParam << _itermin << "\n";
			outfileParam << _dt << "\n";
			outfileParam << _tend << "\n";
			outfileParam << _gravityX << "\n";
			outfileParam << _gravityY << "\n";
			outfileParam << _pr << "\n";
			outfileParam << _beta << "\n";
			outfileParam << _gamma << "\n";
			outfileParam << _q << "\n";
		}
		outfileParam.close();
	}
}

void GeometryGenerator::setFilenames(const char* filenameGeom, const char* filenameParam)
{
	_filenameGeom = filenameGeom;
	_filenameParam = filenameParam;
}

void GeometryGenerator::print() const
{
	std::cout << "==================================================\n";
	std::cout << "GeomGen\n";
	std::cout << "--------------------------------------------------\n";
	std::cout << "Geometry configuration:\n";
	std::cout << "Total size \t=\t" << "(" << _bSizeX-2 << ", " << _bSizeY-2 << ")" << " (without ghost cells)" << "\n";
	std::cout << "Total cells \t=\t" << _totalCells << " (including ghost cells)\n";
	std::cout << "Total length \t=\t(" << _bLengthX << ", " << _bLengthY << ")\n";
	std::cout << "--------------------------------------------------\n";
	std::cout << "Parameter configuration:\n";
	std::cout << "Re \t\t=\t" << _re << std::endl;
	std::cout << "Omega \t\t=\t" << _omega << std::endl;
	std::cout << "Alpha \t\t=\t" << _alpha << std::endl;
	std::cout << "Epsilon \t=\t" << _eps << std::endl;
	std::cout << "Tau \t\t=\t" << _tau << std::endl;
	std::cout << "Itermax \t=\t" << _itermax << std::endl;
	std::cout << "Itermin \t=\t" << _itermin << std::endl;
	std::cout << "dt \t\t=\t" << _dt << std::endl;
	std::cout << "tend \t\t=\t" << _tend << std::endl;
	std::cout << "Gravity \t=\t(" << _gravityX << ", " << _gravityY << ")" << std::endl;
	std::cout << "Pr \t\t=\t" << _pr << std::endl;
	std::cout << "Beta \t\t=\t" << _beta << std::endl;
	std::cout << "Gamma \t\t=\t" << _gamma << std::endl;
	std::cout << "Q \t\t=\t" << _q << std::endl;
	std::cout << "==================================================\n" << std::flush;
}

void GeometryGenerator::initZero()
{
	for (int i=0; i<_bSizeX*_bSizeY; i++){
		_flags[i] = 0;
		_bvu[i] = 0.0;
		_bvv[i] = 0.0;
		_bvp[i] = 0.0;
		_bvt[i] = 0.0;
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
	set(x,y,1,0,0,1,1,0.0,0.0,0.0,0.0);
}

void GeometryGenerator::set(int x, int y, bool obstacle, bool condu, bool condv, bool condp, bool condt, double valu, double valv, double valp, double valt)
{
	char flag = char(obstacle) | (char(condu) << 1) | (char(condv) << 2) | (char(condp) << 3) | (char(condt) << 4);

	set(x,y,flag,valu,valv,valp,valt);
}

void GeometryGenerator::set(int x, int y, char flag, double valu, double valv, double valp, double valt)
{
	int ival = y * (_bSizeX) + x;
	_flags[ival] = flag;
	_bvu[ival] = valu;
	_bvv[ival] = valv;
	_bvp[ival] = valp;
	_bvt[ival] = valt;
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
				double valt = _bvt[ival];
				// right cell
				set(i+1,j,flag,valu,valv,valp,valt);
				// lower cell
				set(i,j-1,flag,valu,valv,valp,valt);
				// lower right cell
				set(i+1,j-1,flag,valu,valv,valp,valt);
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
				_flags[ival] = 1 | 1<<3 | 1<<4; // neumann condition for p and T, dirichlet for u and v
				_bvu[ival] = 1.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = 0.0;
				_bvt[ival] = 0.0;
			} else if (j==0 || i==0 || i==_bSizeX-1){
				// lower, left or right boundary
				_flags[ival] = 1 | 1<<3 | 1<<4; // neumann condition for p and T, dirichlet for u and v
				_bvu[ival] = 0.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = 0.0;
				_bvt[ival] = 0.0;
			}
		}
	}

	setDefaultParameters();

	setFilenames("../data/driven_cavity.geom", "../data/driven_cavity.param");
}

void GeometryGenerator::test_temperature_heating()
{
	// auto balance the number of cells in each dimension
	autoBalance();

	// assure that all values are initialized with zero
	initZero();

	int ival(0);
	for (int i=0; i<_bSizeX; i++){
		for (int j=0; j<_bSizeY; j++){
			ival = j * (_bSizeX) + i;
			if (j==0) {
				// lower boundary
				_flags[ival] = 1 | 1<<3; // neumann condition for p, dirichlet for u and v and T
				_bvu[ival] = 0.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = 0.0;
				_bvt[ival] = 1.0;
			} else if (j==_bSizeY-1 || i==0 || i==_bSizeX-1){
				// upper, left or right boundary
				_flags[ival] = 1 | 1<<3 | 1<<4; // neumann condition for p and T, dirichlet for u and v
				_bvu[ival] = 0.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = 0.0;
				_bvt[ival] = 0.0;
			}
		}
	}

	setDefaultParameters();
	_gravityX = 0.0;
	_gravityY = -9.81;

	setFilenames("../data/test_temperature_heating.geom", "../data/test_temperature_heating.param");
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
				_flags[ival] = 1 | 1<<3 | 1<<4; // neumann condition for p and T, dirichlet for u and v
				_bvu[ival] = 0.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = 0.0;
				_bvt[ival] = 0.0;
			} else if (j==0){
				// lower boundary
				_flags[ival] = 1 | 1<<3 | 1<<4; // neumann condition for p and T, dirichlet for u and v
				_bvu[ival] = 0.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = 0.0;
				_bvt[ival] = 0.0;
			} else if (i==0){
				// left boundary
				_flags[ival] = 1 | 1<<1 | 1<<2 | 1<<4; // dirichlet condition for p, neumann for u and v and T
				_bvu[ival] = 0.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = pressureLeft;
				_bvt[ival] = 0.0;
			} else if (i==_bSizeX-1){
				// right boundary
				_flags[ival] = 1 | 1<<1 | 1<<2 | 1<<4; // dirichlet condition for p, neumann for u and v and T
				_bvu[ival] = 0.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = pressureRight;
				_bvt[ival] = 0.0;
			}
		}
	}

	setDefaultParameters();

	setFilenames("../data/pipe_flow.geom", "../data/pipe_flow.param");
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
			_flags[ival] = 1 | 1<<3 | 1<<4; // neumann condition for p and T, dirichlet for u and v
			_bvu[ival] = 0.0;
			_bvv[ival] = 0.0;
			_bvp[ival] = 0.0;
			_bvt[ival] = 0.0;
		}
	}

	_re = 100.0;

	setFilenames("../data/flow_over_a_step.geom", "../data/flow_over_a_step.param");
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
				_flags[ival] = 1 | 1<<3 | 1<<4; // neumann condition for p and T, dirichlet for u and v
				_bvu[ival] = 0.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = 0.0;
				_bvt[ival] = 0.0;
			}

		}
	}

	fixSingleCells();

	_re = 10000.0;

	setFilenames("../data/karman_vortex_street.geom", "../data/karman_vortex_street.param");
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
				_flags[ival] = 1 | 1<<3 | 1<<4; // neumann condition for p and T, dirichlet for u and v
				_bvu[ival] = 0.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = 0.0;
				_bvt[ival] = 0.0;
			}

		}
	}

	fixSingleCells();

	setDefaultParameters();

	setFilenames("../data/test_case_1.geom", "../data/test_case_1.param");
}

void GeometryGenerator::testCase2()
{
	double xlength = 6.0;
	double ylength = 1.0;
	double alpha = M_PI/4.0;
	double pressureRight = 0.0;
	double pressureLeft = 0.0;
	double width = 0.1;

	karmanVortexStreet(alpha, width, xlength, ylength, pressureLeft, pressureRight);

	// set the left boundary to an inflow boundary
	int ival(0);
	for (int j=0; j<_bSizeY; j++){
		ival = j * (_bSizeX) + 0;
		_flags[ival] = 1 | 1<<3 | 1<<4; // neumann condition for p and T, dirichlet for u and v
		_bvu[ival] = 1.0;
		_bvv[ival] = 0.0;
		_bvp[ival] = 0.0;
		_bvt[ival] = 0.0;
	}

	setDefaultParameters();

	setFilenames("../data/test_case_2.geom", "../data/test_case_2.param");
}

void GeometryGenerator::testCase3()
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
				_flags[ival] = 1 | 1<<4; // neumann condition for p and T, dirichlet for u and v
				_bvu[ival] = 1.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = 0.0;
				_bvt[ival] = 0.0;
			} else if (i==0 || i==_bSizeX-1){
				// left or right boundary
				_flags[ival] = 1 | 1<<4; // neumann condition for p and T, dirichlet for u and v
				_bvu[ival] = 0.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = 0.0;
				_bvt[ival] = 0.0;
			} else if (j==0) {
				_flags[ival] = 1 | 1<<1 | 1<<2 | 1<<4; // neumann condition for p and T, dirichlet for u and v
				_bvu[ival] = 0.0;
				_bvv[ival] = 0.0;
				_bvp[ival] = 0.0;
				_bvt[ival] = 0.0;
			}
		}
	}

	setDefaultParameters();

	setFilenames("../data/test_case_3.geom", "../data/test_case_3.param");
}

void GeometryGenerator::test_twoCellCriterion()
{
	pipeFlow();

	int ival(0);
	for (int i=_bSizeX/2; i<=_bSizeX/2; i++){
		for (int j=0; j<_bSizeY; j++){
			ival = j * (_bSizeX) + i;
			_flags[ival] = 1 | 1<<3 | 1<<4; // neumann condition for p and T, dirichlet for u and v
			_bvu[ival] = 0.0;
			_bvv[ival] = 0.0;
			_bvp[ival] = 0.0;
			_bvt[ival] = 0.0;
		}
	}

	setDefaultParameters();

	setFilenames("../data/test_two_cell_criterion.geom", "../data/test_two_cell_criterion.param");
}

void GeometryGenerator::test_twoCellCriterion2()
{
	test_twoCellCriterion();
	fixSingleCells();
	setFilenames("../data/test_two_cell_criterion_2.geom", "../data/test_two_cell_criterion_2.param");
}

void GeometryGenerator::setDefaultParameters()
{
	_re = 1000.0;
	_omega = 1.7;
	_alpha = 0.9;
	_eps = 0.001;
	_tau = 0.5;
	_itermax = 1000;
	_itermin = 0;
	_dt = 0.5;
	_tend = 30;
	_gravityX = 0.0;
	_gravityY = 0.0;
	_pr = 3.0;
	_beta = 100.0;
	_gamma = 0.9;
	_q = 0.0;
}

void GeometryGenerator::setRe(double re)
{
	_re = re;
}

