#include "CrankNicolson.h"
#include <cmath>
#include <stdexcept>

namespace m2 {

    CrankNicolson::CrankNicolson(const Option& option)
        : option_(option),
        priceMatrix_(option.getTimeDiscr() + 1, option.getSpotDiscr() + 1),
        tridiagonalMatrix_(option.getSpotDiscr() + 1) {
    }

    std::vector<double> CrankNicolson::initializeGrid() const {
        std::vector<double> grid(option_.getSpotDiscr() + 1);
        double dS = 2 * option_.getS0() / option_.getSpotDiscr(); // Grid from 0 to 2 * S0
        for (size_t i = 0; i < grid.size(); ++i) {
            grid[i] = i * dS;
        }
        return grid;
    }

    double CrankNicolson::interpolateRate(double time) const {
        const auto& rates = option_.getRates();
        if (rates.empty()) {
            throw std::runtime_error("No rates provided.");
        }

        for (size_t i = 1; i < rates.size(); ++i) {
            if (time <= rates[i].first) {
                double t1 = rates[i - 1].first, r1 = rates[i - 1].second;
                double t2 = rates[i].first, r2 = rates[i].second;
                return r1 + (time - t1) * (r2 - r1) / (t2 - t1);
            }
        }
        return rates.back().second; // Use the last rate if time exceeds intervals
    }

    void CrankNicolson::applyBoundaryConditions(std::vector<double>& values, double t) const {
        double r = interpolateRate(t);
        double K = option_.getK();
        double S0 = option_.getS0();

        if (option_.getCallPut() == 0) { // Call option
            values[0] = 0.0; // S = 0 -> Payoff is 0
            values.back() = 2 * S0 - K * exp(-r * t); // Far field approximation
        }
        else { // Put option
            values[0] = K * exp(-r * t); // S = 0 -> Payoff is discounted strike
            values.back() = 0.0; // Far field approximation
        }
    }

    void CrankNicolson::solvePDE() {
        auto grid = initializeGrid();
        double dt = option_.getT() / option_.getTimeDiscr();

        // Terminal condition (option payoff at maturity)
        for (size_t i = 0; i < grid.size(); ++i) {
            double payoff = (option_.getCallPut() == 0)
                ? std::max(grid[i] - option_.getK(), 0.0)  // Call payoff
                : std::max(option_.getK() - grid[i], 0.0); // Put payoff
            priceMatrix_(option_.getTimeDiscr(), i) = payoff;
        }

        // Backward time-marching
        for (int t = option_.getTimeDiscr() - 1; t >= 0; --t) {
            double sigma = option_.getSigma();
            double r = interpolateRate(t * dt);
            double dS = grid[1] - grid[0];

            double alpha = 0.5 * dt * (sigma * sigma / (dS * dS) - r / dS);
            double beta = 1.0 + dt * (sigma * sigma / (dS * dS) + r);
            double gamma = -0.5 * dt * (sigma * sigma / (dS * dS) + r / dS);

            // Fill tridiagonal matrix
            for (size_t i = 1; i < grid.size() - 1; ++i) {
                tridiagonalMatrix_.lower(i - 1) = alpha;
                tridiagonalMatrix_.main(i) = beta;
                tridiagonalMatrix_.upper(i - 1) = gamma;
            }

            // Construct RHS (right-hand side)
            std::vector<double> rhs = priceMatrix_.getRow(t + 1);
            applyBoundaryConditions(rhs, t * dt);


            // Solve the tridiagonal system
            std::vector<double> solution = tridiagonalMatrix_.solve(rhs);

            // Store the solution
            for (size_t i = 1; i < grid.size() - 1; ++i) {
                if (option_.getAmerican()) {
                    priceMatrix_(t, i) = std::max(solution[i], priceMatrix_(t + 1, i)); // Early exercise
                }
                else {
                    priceMatrix_(t, i) = solution[i];
                }
            }
        }
    }

    void CrankNicolson::priceOption() {
        solvePDE();
    }

    double CrankNicolson::getOptionPrice() const {
        size_t mid = option_.getSpotDiscr() / 2; // Assuming S0 is the midpoint of the grid
        return priceMatrix_(0, mid);
    }

} // namespace m2
