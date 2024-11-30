#ifndef UTILS_H
#define UTILS_H

#include "Matrix.h"
#include <cmath>
#include <vector>
#include <utility>
#include <stdexcept>

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440 // Value of 1 / sqrt(2)
#endif

// Normal cumulative distribution function
double normalCDF(double x);

// Black-Scholes formula for European call and put options
double blackScholesPrice(bool isCall, double S, double K, double T, double r, double sigma);

// Maximum between two numbers
float max(float a, float b);

// Linear interpolation of interest rate
float interpolateRate(float t, const std::vector<std::pair<float, float>>& rates);

// Overload the += operator for std::vector<float>
std::vector<float>& operator+=(std::vector<float>& lhs, const std::vector<float>& rhs);

// Using the Crout algorithm to calculate Upper and Lower decomposition of T2 and solve the system
void crout(m2::Matrix& T2, std::vector<float>& W, std::vector<float>& V, m2::Matrix& result, unsigned int M);

// Compute the average rate over time using the piecewise rate structure
double computeAverageRate(const std::vector<std::pair<double, double>>& rates, double T);

#endif // !UTILS_H
