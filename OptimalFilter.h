#ifndef _OptimalFilter_OptimalFilter_h
#define _OptimalFilter_OptimalFilter_h

#include <Core/Core.h>
#include <CtrlLib/CtrlLib.h>
#include <plugin/Eigen/Eigen.h>
#include <unsupported/Eigen/Polynomials>
#include <ScatterCtrl/ScatterCtrl.h>
#include <complex>

using namespace Upp;
using namespace Eigen;
using namespace std;

#define LAYOUTFILE <OptimalFilter/OptimalFilter.lay>
#include <CtrlCore/lay.h>
#include <Draw/Draw.h>

// min/max
#define MAX2(xi_a, xi_b) ((xi_a >= xi_b) ? (xi_a) : (xi_b))
#define MIN2(xi_a, xi_b) ((xi_a >= xi_b) ? (xi_b) : (xi_a))

// math constants
static const double pi = 3.1415926535897932384626433832795,
				    r2d = 57.295779513082320876798154814105;
static const std::complex<double> I(0.0, 1.0);

// circle 
Pointf circleDraw(double t) {return Pointf(cos(2*pi*t), sin(2*pi*t));}

class OptimalFilter : public WithOptimalFilterLayout<TopWindow> {
public:

	// constructor
	typedef OptimalFilter CLASSNAME;
	OptimalFilter();
		
	// filter order
	void filterOrderSpinerFunc();
	
	// button callbacks
	void calculateBtnFunc();
	void aboutBtnFunc();
	
	// optimal filter calculation
	void OptimalFilterCalc(const int xi_order, const double xi_fc, const double);

	
private:
	VectorMap<double, double> factStable;
	VectorMap<double, double> zero;
	double xMin, xMax, yMin, yMax;
	
	Vector<double> magnitude, phase, Fs;
	double magMin, magMax, phaseMin, phaseMax;
	//VectorMap<double, double> magCross;
	
	VectorMap<double, double> step;//, impulse;
	double stepMax, stepMin;//, impulseMin, impulseMax;
};

#endif
