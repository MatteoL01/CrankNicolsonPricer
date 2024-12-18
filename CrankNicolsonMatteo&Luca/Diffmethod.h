#ifndef DIFFMETHOD_H
#define DIFFMETHOD_H
#include "Matrix.h"
#include "Option.h"
#include "Utils.h"
#include <iostream>
#include <cmath>
#include <vector>

namespace m2
{

	/**
	 * @class American
	 * @brief Represents an American option and provides methods for pricing and calculating Greeks.
	 */
	class American : public Option
	{
	private:
		Matrix values_; ///< Matrix storing the option values.
		double dt_; ///< Time step size.
		double ds_; ///< Spot price step size.
		double Smax_; ///< Maximum spot price for calculations.

		double price_; ///< Final price of the American option.
		std::vector<double> T0prices_; ///< Option prices at initial time T0 for varying spot prices.
		double optionPrice_; ///< Option price.

		std::vector<double> delta_; ///< Delta values for the option.
		std::vector<double> gamma_; ///< Gamma values for the option.
		std::vector<double> theta_; ///< Theta values for the option.
		double vega_; ///< Vega value for the option.
		double rho_; ///< Rho value for the option.

		std::vector<double> boundary_; ///< Boundary conditions for the option.
	public:
		/**
		 * @brief Constructs an American option using the given Option object.
		 * @param opt Option object with parameters.
		 */
		American(Option & opt);
		
		/**
		 * @brief Prices the American call option.
		 */
		void priceCall();
		/**
		 * @brief Prices the American put option.
		 */
		void pricePut();

		/**
		 * @brief Calculates the Delta (rate of change of price with respect to underlying).
		 */
		void calculateDelta();

		/**
		 * @brief Calculates the Gamma (rate of change of Delta with respect to underlying).
		 */
		void calculateGamma();

		/**
		 * @brief Calculates the Theta (rate of change of price with respect to time).
		 */
		void calculateTheta();

		/**
		 * @brief Calculates the Vega (rate of change of price with respect to volatility).
		 */
		void calculateVega();

		/**
		 * @brief Calculates the Rho (rate of change of price with respect to interest rate).
		 */
		void calculateRho();

		// Getters

		/**
		 * @brief Retrieves the matrix of option values.
		 * @return Matrix of values.
		 */
		Matrix get_val() const { return values_; }

		/**
		 * @brief Retrieves the final price of the American option.
		 * @return Price of the option.
		 */
		double getPrice() const { return price_;}

		/**
		 * @brief Retrieves option prices at T0 for varying spot prices.
		 * @return Vector of prices.
		 */
		std::vector<double> getT0prices() const { return T0prices_; }

		// Greeks
		std::vector<double> getDeltaGraph() const { return delta_; }
		double getDelta() const { return delta_[S0_ / ds_]; }
		std::vector<double> getGammaGraph() const { return gamma_; }
		double getGamma() const { return gamma_[S0_ / ds_]; }
		std::vector<double> getThetaGraph() const { return theta_; }
		double getTheta() const { return theta_[S0_ / ds_]; }
		double getVega() const { return vega_; }
		double getRho() const { return rho_; }

		// Boundaries
		std::vector<double> getBound() const { return boundary_; }

		/**
		 * @brief Prints the matrix of option values to the console.
		 */
		void printMatrix() const {
			std::cout << "Matrix Values: \n" << values_ << std::endl;
		}

	};

	/**
	 * @class European
	 * @brief Represents a European option and provides methods for pricing and calculating Greeks.
	 */
	class European : public Option
	{
	private:
		Matrix values_; ///< Matrix storing the option values.
		double dt_; ///< Time step size.
		double ds_; ///< Spot price step size.
		double Smax_; ///< Maximum spot price for calculations.

		double price_; ///< Final price of the European option.
		std::vector<double> T0prices_; ///< Option prices at initial time T0 for varying spot prices.
		double optionPrice_; ///< Option price.

		std::vector<double> delta_; ///< Delta values for the option.
		std::vector<double> gamma_; ///< Gamma values for the option.
		std::vector<double> theta_; ///< Theta values for the option.
		double vega_; ///< Vega value for the option.
		double rho_; ///< Rho value for the option.

		std::vector<double> boundary_; ///< Boundary conditions for the option.
	public:
		/**
		 * @brief Constructs a European option using the given Option object.
		 * @param opt Option object with parameters.
		 */
		European(Option& opt);

		/**
		 * @brief Prices the European call option.
		 */
		void priceCall();
		/**
		 * @brief Prices the European put option.
		 */
		void pricePut();

		/**
		 * @brief Calculates the Delta (rate of change of price with respect to underlying).
		 */
		void calculateDelta();

		/**
		 * @brief Calculates the Gamma (rate of change of Delta with respect to underlying).
		 */
		void calculateGamma();

		/**
		 * @brief Calculates the Theta (rate of change of price with respect to time).
		 */
		void calculateTheta();

		/**
		 * @brief Calculates the Vega (rate of change of price with respect to volatility).
		 */
		void calculateVega();

		/**
		 * @brief Calculates the Rho (rate of change of price with respect to interest rate).
		 */
		void calculateRho();

		// Getters
		/**
		 * @brief Retrieves the matrix of option values.
		 * @return Matrix of values.
		 */
		Matrix get_val() const { return values_; }

		/**
		 * @brief Retrieves the final price of the American option.
		 * @return Price of the option.
		 */
		double getPrice() const { return price_; }

		/**
		 * @brief Retrieves option prices at T0 for varying spot prices.
		 * @return Vector of prices.
		 */
		std::vector<double> getT0prices() const { return T0prices_; }

		// Greeks
		std::vector<double> getDeltaGraph() const { return delta_; }
		double getDelta() const { return delta_[S0_ / ds_]; }
		std::vector<double> getGammaGraph() const { return gamma_; }
		double getGamma() const { return gamma_[S0_ / ds_]; }
		std::vector<double> getThetaGraph() const { return theta_; }
		double getTheta() const { return theta_[S0_ / ds_]; }
		double getVega() const { return vega_; }
		double getRho() const { return rho_; }

		// Boundaries
		std::vector<double> getBound() const { return boundary_; }

		/**
		 * @brief Prints the matrix of option values to the console.
		 */
		void printMatrix() const {
			std::cout << "Matrix Values: \n" << values_ << std::endl;
		}
	};
}

#endif // !DIFFMETHOD_H

