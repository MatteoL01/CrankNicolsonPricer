#ifndef OPTION_H
#define OPTION_H

#include <vector>
#include <utility>
#include <string>


namespace m2 {
    class Option {
    protected:
        bool call_; // 1 for call, 0 for put
        bool american_; // 0 for European, 1 for American
        double T_;           // Maturity date
        double K_;           // Strike price
        unsigned int timeDiscr_; // Discretization in time
        unsigned int spotDiscr_; // Discretization in spot
        double S0_;         // Initial price of the underlying
        double sigma_;       // Volatility
        unsigned int time_steps_; // Interest rate steps
        std::vector<std::pair<double, double>> rates_;  // Interest rate as (time, rate)

        void validateInput(); // Validate parameters

    public:
        // Default constructor
        Option();

        // Open input helper
        void openInput(const std::string& filePath);

        // Debug/Print function
        void printParameters() const;

        //double getInterestRate(const std::vector<std::pair<double, double>>& rates, double time) const;

        // Getter Functions
        bool getCallPut() const { return call_; }
        bool getAmerican() const { return american_; }
        double getT() const { return T_; }
        double getK() const { return K_; }
        double getS0() const { return S0_; }
        double getSigma() const { return sigma_; }
        unsigned int getTimeDiscr() const { return timeDiscr_; }
        unsigned int getSpotDiscr() const { return spotDiscr_; }
        unsigned int getTimeSteps() const { return time_steps_; }
        const std::vector<std::pair<double, double>>& getRates() const { return rates_; }
    };

    
}

#endif // OPTION_H
