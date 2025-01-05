/**
 * @file utils.h
 * @brief Utility functions and constants for financial option pricing and mathematical operations.
 *
 * This header provides various mathematical utilities, constants, and helper functions
 * for implementing financial models and numerical methods.
 */

#ifndef UTILS_H
#define UTILS_H

#include "Diffmethod.h"
#include "Matrix.h"
#include <cmath>
#include <vector>
#include <utility>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

 /** @def M_PI
  *  @brief Defines the value of pi if not already defined.
  */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

  /** @def M_SQRT1_2
   *  @brief Defines the value of 1/sqrt(2) if not already defined.
   */
#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

namespace m2 {

	/**
	 * @brief Computes the cumulative distribution function of the standard normal distribution.
	 * @param x The input value.
	 * @return The probability that a standard normal random variable is less than or equal to x.
	 */
	double normalCDF(double x);

	/**
	 * @brief Computes the probability density function of the standard normal distribution.
	 * @param d The input value.
	 * @return The value of the standard normal density function at d.
	 */
	double normalPDF(double d);

	/**
	 * @brief Computes the price of a European option using the Black-Scholes formula.
	 * @param isCall True for a call option, false for a put option.
	 * @param S The current price of the underlying asset.
	 * @param K The strike price of the option.
	 * @param T The time to maturity (in years).
	 * @param r The risk-free interest rate.
	 * @param sigma The volatility of the underlying asset.
	 * @return The price of the option.
	 */
	double blackScholesPrice(bool isCall, double S, double K, double T, double r, double sigma);

	/**
	 * @brief Computes the maximum of two floating-point numbers.
	 * @param a The first number.
	 * @param b The second number.
	 * @return The larger of the two numbers.
	 */
	double max(double a, double b);

	/**
	 * @brief Linearly interpolates the interest rate at a given time.
	 * @param t The time for which to interpolate the rate.
	 * @param rates A vector of (time, rate) pairs representing the piecewise interest rate structure.
	 * @return The interpolated rate.
	 */
	double interpolateRate(double t, const std::vector<std::pair<double, double>>& rates);

	/**
	 * @brief Adds the elements of one vector to another (element-wise).
	 * @param lhs The left-hand-side vector to which elements are added.
	 * @param rhs The right-hand-side vector from which elements are added.
	 * @return A reference to the updated lhs vector.
	 */
	std::vector<double>& operator+=(std::vector<double>& lhs, const std::vector<double>& rhs);

	/**
	 * @brief Solves a tridiagonal system of equations using the Crout decomposition.
	 * @param T2 The tridiagonal matrix representing the system.
	 * @param W The right-hand-side vector of the system.
	 * @param V The solution vector (updated in-place).
	 * @param M The size of the system (number of equations).
	 */
	void crout(const Matrix& T2, const std::vector<double>& W, std::vector<double>& V, long int M);

	/**
	 * @brief Computes the average interest rate over a given time period using a piecewise rate structure.
	 * @param rates A vector of (time, rate) pairs representing the piecewise interest rate structure.
	 * @param T The time period over which to compute the average rate.
	 * @return The average rate.
	 */
	double computeAverageRate(const std::vector<std::pair<double, double>>& rates, double T);

	/**
	 * @brief Writes the results of the option pricing calculations to a text file.
	 * @param price The calculated option price.
	 * @param delta The option delta.
	 * @param gamma The option gamma.
	 * @param theta The option theta.
	 * @param vega The option vega.
	 * @param rho The option rho.
	 * @param T0prices A vector of option prices at time 0.
	 * @param T0deltas A vector of option deltas at time 0.
	 * @param boundaries A vector of boundary conditions used in the calculations.
	 */
	void writeOutputTxt(double price, double delta, double gamma, double theta, double vega, double rho, const std::vector<double>& T0prices, const std::vector<double>& T0deltas, const std::vector<double>& boundaries);

} // namespace m2

#endif // UTILS_H
