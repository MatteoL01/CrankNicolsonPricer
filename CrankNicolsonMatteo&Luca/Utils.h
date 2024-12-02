#ifndef UTILS_H
#define UTILS_H
#include "Diffmethod.h"
#include "Matrix.h"
#include <cmath>
#include <vector>
#include <utility>
#include <stdexcept>

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440 // Value of 1 / sqrt(2)
#endif

namespace m2 {
	// Normal cumulative distribution function
	double normalCDF(double x);

	// Black-Scholes formula for European call and put options
	double blackScholesPrice(bool isCall, double S, double K, double T, double r, double sigma);

	// Maximum between two numbers
	float max(float a, float b);

	// Linear interpolation of interest rate
	float interpolateRate(double t, const std::vector<std::pair<double, double>>& rates);

	// Overload the += operator for std::vector<float>
	std::vector<double>& operator+=(std::vector<double>& lhs, const std::vector<double>& rhs);

	// Using the Crout algorithm to calculate Upper and Lower decomposition of T2 and solve the system
	void crout(Matrix& T2, std::vector<double>& W, std::vector<double>& V, int M);

	// Compute the average rate over time using the piecewise rate structure
	double computeAverageRate(const std::vector<std::pair<double, double>>& rates, double T);

	// Write the results of the calculations in a txt file 
	void writeOutputTxt();
}

#endif // !UTILS_H

