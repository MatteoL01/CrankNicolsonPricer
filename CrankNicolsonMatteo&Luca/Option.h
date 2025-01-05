#ifndef OPTION_H
#define OPTION_H

#include <vector>
#include <utility>
#include <string>

namespace m2 {

    /**
     * @class Option
     * @brief Represents a financial option and provides utilities for option pricing.
     */
    class Option {
    protected:
        bool call_; ///< True for a call option, false for a put option.
        bool american_; ///< True for an American option, false for a European option.
        double T_; ///< Maturity date of the option (in years).
        double K_; ///< Strike price of the option.
        long int N_; ///< Number of time discretization steps.
        long int M_; ///< Number of spot discretization steps.
        double S0_; ///< Initial price of the underlying asset.
        double sigma_; ///< Volatility of the underlying asset.
        long int time_steps_; ///< Number of steps for interest rate discretization.
        std::vector<std::pair<double, double>> rates_; ///< Interest rates as pairs of (time, rate).

        /**
         * @brief Validates the input parameters for the option.
         * @throw std::invalid_argument if parameters are invalid.
         */
        void validateInput();

    public:
        /**
         * @brief Default constructor.
         * Initializes the option with default parameters.
         */
        Option();

        /**
         * @brief Loads option parameters from an input file.
         * @param filePath Path to the input file containing option parameters.
         */
        void openInput(const std::string& filePath);

        /**
         * @brief Prints the option parameters to the standard output.
         */
        void printParameters() const;

        // Getter Functions

        /**
         * @brief Checks whether the option is a call or put.
         * @return True for a call option, false for a put option.
         */
        bool getCallPut() const { return call_; }


        /**
         * @brief Checks whether the option is American or European.
         * @return True for an American option, false for a European option.
         */
        bool getAmerican() const { return american_; }

        /**
         * @brief Gets the maturity date of the option.
         * @return Maturity date in years.
         */
        double getT() const { return T_; }

        /**
         * @brief Gets the strike price of the option.
         * @return Strike price.
         */
        double getK() const { return K_; }

        /**
         * @brief Gets the initial price of the underlying asset.
         * @return Initial price of the underlying.
         */
        double getS0() const { return S0_; }

        /**
         * @brief Gets the volatility of the underlying asset.
         * @return Volatility as a double.
         */
        double getSigma() const { return sigma_; }

        /**
         * @brief Gets the number of time discretization steps.
         * @return Number of time steps.
         */
        long int getTimeDiscr() const { return N_; }

        /**
         * @brief Gets the number of spot discretization steps.
         * @return Number of spot steps.
         */
        long int getSpotDiscr() const { return M_; }

        /**
         * @brief Gets the number of time steps for interest rate discretization.
         * @return Number of interest rate steps.
         */
        long int getTimeSteps() const { return time_steps_; }

        /**
         * @brief Gets the interest rate data.
         * @return A vector of pairs, where each pair contains (time, rate).
         */
        const std::vector<std::pair<double, double>>& getRates() const { return rates_; }
    };

}

#endif // OPTION_H
