#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
using namespace std;

// this class contains the parameters for the genes in the population, such
// as permissible ranges of each gene
class Genetics
{
public:
	vector<double> mins;		// minimum value a gene can posess
	vector<double> maxes;	// maximum value a gene can posess
	vector<double> range;	// maxes - mins

	vector< vector<double> > targets;	// the target data set
};

// this class contains information about a specific indivudual, such as the
// value for each gene and the performance of the individual
class Individual
{
public:
	// the gene values for this individual
	vector<double> gene;

	// the sensitivity of each gene; how much change in performance results from 
	// a small (1% of range) change in the corresponding gene value at this 
	// individual's value
	vector<double> sensitivity;

	// the overall performance of this individual
	double performance;
	double stability;
	
};

// this class owns the population, and handles the evolution
class Population
{
public:
	vector<Individual> population;
	Genetics genetics;
	int generation;

	// Generate a new population with random genetic values
	void InitializePopulation(int members);

	// For each generation:

	// run simulations on all members
	// sort individuals by performance
	// preserve the genes of the top performer(s) for breeding
	// To create next generation:
	// create some all-new individuals
	// mutate top performers of this generation
	// mix genes of this generation's top perfomers with genes of overall top performers
};

// Cubic Bezier curve, generating a curve from (0, 0) to (1, 1)
// (a, b) is the lower control point, (b, c) is the upper control point
// t is the parameter, and should range from [0 - 1].
void BezierCurve(double a, double b, double c, double d, double t, double &xout, double &yout)
{
	// lower left endpoint
	double P1x = 0;
	double P1y = 0;
	// lower control endpoint
	double P2x = a;
	double P2y = b;
	// upper control endpoint
	double P3x = c;
	double P3y = d;
	// upper right endpoint
	double P4x = 1;
	double P4y = 1;

	xout = P1x * pow(1 - t, 3) + P2x * 3 * t * pow(1 - t, 2) + 
		P3x * 3 * pow(t, 2) * (1 - t) + P4x * pow(t, 3);
	yout = P1y * pow(1 - t, 3) + P2y * 3 * t * pow(1 - t, 2) + 
		P3y * 3 * pow(t, 2) * (1 - t) + P4y * pow(t, 3);

	return;
}


// generates a 256 element array with a cumulative amount of powder burned from 0 to 1
// x and y are the coordinates of the lower control point, z is the magnitude of the 
// upper control point (which is always horiztonal, to ensyre the curve approaches
// 1.0 tangientally
void BurnCurve(double x, double y, double z, vector<double> &array)
{
	int i = -1;
	double t;
	double last;

	// set up array and initialize it to invalid values
	array.resize(256);
	for(i = 0; i < 256; ++i)
	{
		array[i] = -1;
	}

	// since this is a parametric curve, sweep the parameter and get the X and Y values
	for(t = 0; t < 1; t += .01)
	{
		double dx, dy;

		// set Y of second control point to 1, so the curve approaches 1 tangientally
		BezierCurve(x, y,  1 - z, 1, t, dx, dy);
		
		// clip values to the range 0-1
		if(dx > 1) dx = 1;
		if(dy > 1) dy = 1;
		if(dx < 0) dx = 0;
		if(dy < 0) dy = 0;

		int index = dx * 255;

		// if there is not already a value at this point, enter it
		if(array[(long)index] == -1)	array[(long)index] = dy;
	}

	// the array may be sparsely filled, so fill in blank spots
	// make sure no value is less than the previous value
	last = 0;
	for(i = 0; i < 256; ++i)
	{
		if(array[i] < last)
			array[i] = last;
		else
			last = array[i];
	}
	return;
}

// parameters: 
// 0 - 2:	a, b, c parameters of burn curve, range 0 to 1
// 3:		time for full burn in ms at 10kpsi, values around 1
// 4:		change in burn time for pressure changes, values around 1 
// 5:		cubic inches of gas produced per grain of propellant at 1 psi, values around 300

/* various constants
 * density of nitrocellulose, 1.23 g/ml  http://www.chemicalbook.com/ChemicalProductProperty_EN_CB6781086.htm
 * Isenotropic expansion of a gas, from https://ccrma.stanford.edu/~jos/pasp/Adiabatic_Gas_Constant.html
 * pressure = pressure * (volume 1 / volume 2) ^ Y, where Y is 1.6 for monatomic (noble gases), 1.4 for diatomic (N2), 
 * and 1.28 for triatomic (CO2).  A value of about 1.35 would probably work for nitrocellulose combustion gases.
*/

// takes the input parameters, runs the simulation, returns the outputs
void RunSimulation(vector<double> &parameters, vector<double> &outputs)
{
	// generate the burn curve for this run
	vector<double> curve;
	BurnCurve(parameters[0], parameters[1], parameters[2], curve);
	
}

int main(int argc, char **argv)
{
	vector<double> curve;
	BurnCurve(0.5, 0.0, 0.5, curve);

	for(int y = 0; y < curve.size(); y += 8)
	{
		printf("%f %f ", (double)y / 256.0, curve[y]);
		for(int x = 0; x < curve[y] * 100; ++x)
		{
			printf("*");
		}
		printf("\n");
	}	

	return 0;
}


/*
 *  example output of a burn curve:
 
	vector<double> curve;
	BurnCurve(0.5, 0.0, 0.5, curve);

	for(int y = 0; y < curve.size(); y += 8)
	{
		printf("%f %f ", (double)y / 256.0, curve[y]);
		for(int x = 0; x < curve[y] * 100; ++x)
		{
			printf("*");
		}
		printf("\n");
	}	

0.000000 0.000000 
0.031250 0.001184 *
0.062500 0.004672 *
0.093750 0.014014 **
0.125000 0.022842 ***
0.156250 0.039744 ****
0.187500 0.053312 ******
0.218750 0.076874 ********
0.250000 0.104000 ***********
0.281250 0.134366 **************
0.312500 0.179334 ******************
0.343750 0.216000 **********************
0.375000 0.268192 ***************************
0.406250 0.323456 *********************************
0.437500 0.381024 ***************************************
0.468750 0.440128 *********************************************
0.500000 0.500000 ***************************************************
0.531250 0.559872 ********************************************************
0.562500 0.633542 ****************************************************************
0.593750 0.690606 **********************************************************************
0.625000 0.731808 **************************************************************************
0.656250 0.784000 *******************************************************************************
0.687500 0.832352 ************************************************************************************
0.718750 0.865634 ***************************************************************************************
0.750000 0.896000 ******************************************************************************************
0.781250 0.923126 *********************************************************************************************
0.812500 0.946688 ***********************************************************************************************
0.843750 0.960256 *************************************************************************************************
0.875000 0.977158 **************************************************************************************************
0.906250 0.985986 ***************************************************************************************************
0.937500 0.995328 ****************************************************************************************************
0.968750 0.998816 ****************************************************************************************************

*/
