#include "Option.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

namespace m2 {

    // Default constructor
    Option::Option()
        : call_(0), american_(0), T_(0.0), K_(0.0), N_(0), M_(0),
        S0_(0.0), sigma_(0.0), time_steps_(0) {
    }


    // Validate input parameters
    void Option::validateInput() {
        if (T_ <= 0) throw std::invalid_argument("Maturity (T) must be greater than 0.");
        if (K_ <= 0) throw std::invalid_argument("Strike price (K) must be greater than 0.");
        if (sigma_ <= 0) throw std::invalid_argument("Volatility (sigma) must be greater than 0.");
        if (N_ == 0) throw std::invalid_argument("Time discretization must be greater than 0.");
        if (M_ == 0) throw std::invalid_argument("Spot discretization must be greater than 0.");
        if (rates_.empty()) throw std::invalid_argument("Rates vector cannot be empty.");
    }

    // Open input implementation
    void Option::openInput(const std::string& filePath) {
        std::ifstream file(filePath);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open file: " + filePath);
        }

        try {
            std::string line;

            // Read contract type
            std::getline(file, line);
            if (line == "call") call_ = true;
            else if (line == "put") call_ = false;
            else throw std::invalid_argument("The contract type should be 'call' or 'put'.");

            // Read exercise type
            std::getline(file, line);
            if (line == "american") american_ = true;
            else if (line == "european") american_ = false;
            else throw std::invalid_argument("The exercise type should be 'american' or 'european'.");

            // Read other parameters
            std::getline(file, line); T_ = std::stod(line);
            std::getline(file, line); K_ = std::stod(line);
            std::getline(file, line); N_ = std::stoi(line);
            std::getline(file, line); M_ = std::stoi(line);
            std::getline(file, line); S0_ = std::stod(line);
            std::getline(file, line); sigma_ = std::stod(line);
            std::getline(file, line); time_steps_ = std::stoul(line);

            // Clear and read rates
            rates_.clear();
            while (std::getline(file, line)) {
                std::istringstream stream(line);
                double time, rate;
                if (stream >> time >> rate) {
                    rates_.emplace_back(time, rate);
                }
                else {
                    throw std::invalid_argument("Invalid rate pair format in input file.");
                }
            }

            file.close();
            validateInput(); // Validate all parameters after reading
        }
        catch (const std::exception& e) {
            file.close();
            throw std::runtime_error("Error parsing input file: " + std::string(e.what()));
        }
    }

    // Debug/Print function
    void Option::printParameters() const {
        std::cout << "Contract Type: " << (call_ ? "Call" : "Put") << "\n";
        std::cout << "Exercise Type: " << (american_ ? "American" : "European") << "\n";
        std::cout << "Maturity (T): " << T_ << "\n";
        std::cout << "Strike Price (K): " << K_ << "\n";
        std::cout << "Time Discretization: " << N_ << "\n";
        std::cout << "Spot Discretization: " << M_ << "\n";
        std::cout << "Initial Price (S0): " << S0_ << "\n";
        std::cout << "Volatility (sigma): " << sigma_ << "\n";
        std::cout << "Time Steps: " << time_steps_ << "\n";
        std::cout << "Rates:\n";
        for (const auto& rate : rates_) {
            std::cout << "  Time: " << rate.first << ", Rate: " << rate.second << "\n";
        }
    }

    
}