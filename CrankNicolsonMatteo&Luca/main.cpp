#include "Option.h"
#include "CrankNicolson.h"
#include "utils.h"  
#include "Matrix.h"
#include <iostream>
#include <iomanip> // For formatting

int main() {
    try {
        // Load option parameters
        m2::Option opt;
        opt.openInput("C:/Users/matte/Documents/option_data.txt");

        // Create Crank-Nicholson pricer
        m2::CrankNicolson pricer(opt);

        // Price the option
        pricer.priceOption();

        m2::Matrix priceMatrix = pricer.getPriceMatrix();

        priceMatrix.print();

        // Get the option price from Crank-Nicholson
        double crankNicholsonPrice = pricer.getOptionPrice();
        std::cout << "Crank-Nicholson Price: " << crankNicholsonPrice << "\n";

        // Compute average interest rate over the period
        double T = opt.getT();
        const auto& rates = opt.getRates();  // Get rates from the Option object
        double avgRate = computeAverageRate(rates, T);

        std::cout << "Average Rate: " << avgRate << "\n";

        // Get the price from Black-Scholes formula
        bool isCall = (opt.getCallPut() == 0);  // 0 for call, 1 for put
        double S = opt.getS0();
        double K = opt.getK();
        double sigma = opt.getSigma();

        double blackScholesPriceValue = blackScholesPrice(isCall, S, K, T, avgRate, sigma);
        std::cout << "Black-Scholes Price: " << blackScholesPriceValue << "\n";

        // Compare the results
        double error = std::abs(crankNicholsonPrice - blackScholesPriceValue);
        std::cout << "Difference: " << error << "\n";

    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
