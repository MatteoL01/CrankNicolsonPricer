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
	class American
	{
	private:
		Matrix values_;
		bool call_; // 1 for call, 0 for put
		double T_;           // Maturity date
		double K_;           // Strike price
		unsigned int N_; // Discretization in time
		unsigned int M_; // Discretization in spot
		double S0_;         // Initial price of the underlying
		double sigma_;       // Volatility
		std::vector<std::pair<double, double>> rates_;  // Interest rate as (time, rate)

		double price_;

	public:
		American(Option & opt);
		
		void priceCall();
		void pricePut();

		Matrix get_val() { return values_; }

		double getPrice() const { return price_;}

		void printMatrix() {
			std::cout << "Matrix Values: \n" << values_ << std::endl;
		}

	};

	class European
	{
	private:
		Matrix values_;
		bool call_; // 1 for call, 0 for put
		double T_;           // Maturity date
		double K_;           // Strike price
		unsigned int N_; // Discretization in time
		unsigned int M_; // Discretization in spot
		double S0_;         // Initial price of the underlying
		double sigma_;       // Volatility
		std::vector<std::pair<double, double>> rates_;  // Interest rate as (time, rate)

		double price_;
	public:
		European(Option& opt);

		void priceCall();
		void pricePut();

		Matrix get_val() {return values_; }

		double getPrice() const { return price_; }

		void printMatrix() {
			std::cout << "Matrix Values: \n" << values_ << std::endl;
		}
	};
}

#endif // !DIFFMETHOD_H

