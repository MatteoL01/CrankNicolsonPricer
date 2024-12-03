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

		// store the values of the Greeks
		std::vector<double> delta_;
	public:
		American(Option & opt);
		
		void priceCall();
		void pricePut();

		void calculateDelta();

		Matrix get_val() const { return values_; }

		double getPrice() const { return price_;}

		std::vector<double> getT0prices() const { return T0prices_; }

		// greeks
		std::vector<double> getDeltaGraph() const { return delta_; }
		double getDelta() const { return delta_[S0_ / ds_]; }

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

		// store the values of the Greeks
		std::vector<double> delta_;
	public:
		European(Option& opt);

		void priceCall();
		void pricePut();

		void calculateDelta();

		Matrix get_val() const {return values_; }

		double getPrice() const { return price_; }

		std::vector<double> getT0prices() const { return T0prices_; }

		// greeks
		std::vector<double> getDeltaGraph() const { return delta_; }
		double getDelta() const { return delta_[S0_ / ds_]; }

		void printMatrix() const {
			std::cout << "Matrix Values: \n" << values_ << std::endl;
		}
	};
}

#endif // !DIFFMETHOD_H

