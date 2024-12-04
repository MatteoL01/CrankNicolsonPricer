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
	class American : public Option
	{
	private:
		Matrix values_;
		double dt_;
		double ds_;
		double Smax_;
		
		// store the values of the prices
		double price_; // final price of American Option
		std::vector<double> T0prices_; // price of Option at time T0 varying S
		double optionPrice_;

		// store the values of the Greeks
		std::vector<double> delta_;
		std::vector<double> gamma_;
		std::vector<double> theta_;
		double vega_;
		double rho_;

		// Store the values of Boundary conditions
		std::vector<double> boundary_;
	public:
		American(Option & opt);
		
		void priceCall();
		void pricePut();

		void calculateDelta();
		void calculateGamma();
		void calculateTheta();
		void calculateVega();
		void calculateRho();

		Matrix get_val() const { return values_; }

		double getPrice() const { return price_;}

		std::vector<double> getT0prices() const { return T0prices_; }

		// greeks
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

		void printMatrix() const {
			std::cout << "Matrix Values: \n" << values_ << std::endl;
		}

	};

	class European : public Option
	{
	private:
		Matrix values_;
		double dt_;
		double ds_;
		double Smax_;

		// store the values of the prices
		double price_; // final price of European Option
		std::vector<double> T0prices_; // price of Option at time T0 varying S
		double optionPrice_;

		// store the values of the Greeks
		std::vector<double> delta_;
		std::vector<double> gamma_;
		std::vector<double> theta_;
		double vega_;
		double rho_;

		// Store the values of Boundary conditions
		std::vector<double> boundary_;
	public:
		European(Option& opt);

		void priceCall();
		void pricePut();

		void calculateDelta();
		void calculateGamma();
		void calculateTheta();
		void calculateVega();
		void calculateRho();

		Matrix get_val() const {return values_; }

		double getPrice() const { return price_; }

		std::vector<double> getT0prices() const { return T0prices_; }

		// greeks
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

		void printMatrix() const {
			std::cout << "Matrix Values: \n" << values_ << std::endl;
		}
	};
}

#endif // !DIFFMETHOD_H

