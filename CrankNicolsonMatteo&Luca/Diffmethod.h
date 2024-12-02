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
		/*
		bool call_; // true for call, 0 for put
		double T_;           // Maturity date
		double K_;           // Strike price
		unsigned int N_; // Discretization in time
		unsigned int M_; // Discretization in spot
		double S0_;         // Initial price of the underlying
		double sigma_;       // Volatility
		std::vector<std::pair<double, double>> rates_;  // Interest rate as (time, rate)
		*/

		// store the values of the prices
		double price_; // final price of American Option
		std::vector<double> T0prices_; // price of Option at time T0 varying S

		// store the values of the Greeks
		Matrix delta_;
	public:
		American(Option & opt);
		
		void priceCall();
		void pricePut();

		void calculateDelta();

		Matrix get_val() const { return values_; }

		double getPrice() const { return price_;}

		std::vector<double> getT0prices() const { return T0prices_; }

		void printMatrix() const {
			std::cout << "Matrix Values: \n" << values_ << std::endl;
		}

	};

	class European : public Option
	{
	private:
		Matrix values_;
		
		// store the values of the prices
		double price_; // final price of European Option
		std::vector<double> T0prices_; // price of Option at time T0 varying S

		// store the values of the Greeks
		Matrix delta_;
	public:
		European(Option& opt);

		void priceCall();
		void pricePut();

		void calculateDelta();

		Matrix get_val() const {return values_; }

		double getPrice() const { return price_; }

		std::vector<double> getT0prices() const { return T0prices_; }

		void printMatrix() const {
			std::cout << "Matrix Values: \n" << values_ << std::endl;
		}
	};
}

#endif // !DIFFMETHOD_H

