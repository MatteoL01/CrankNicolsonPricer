#ifndef CRANKNICOLSON_H
#define CRANKNICOLSON_H

#include "Option.h"
#include "Matrix.h"
#include "TridiagonalMatrix.h"
#include <vector>

namespace m2 {

    class CrankNicolson {
    private:
        const Option& option_;                  // Reference to the Option object
        Matrix priceMatrix_;                    // Stores option prices at each grid point
        TridiagonalMatrix tridiagonalMatrix_;   // Tridiagonal matrix for solving the system

        // Helper functions
        std::vector<double> initializeGrid() const;  // Initializes the spot price grid
        double interpolateRate(double time) const;   // Interpolates the risk-free rate

        // Boundary condition handlers
        void applyBoundaryConditions(std::vector<double>& values, double t) const;

        // PDE solver
        void solvePDE(); // Implements the Crank-Nicholson scheme

    public:
        explicit CrankNicolson(const Option& option);

        // Main pricing method
        void priceOption();

        // Accessors
        const Matrix& getPriceMatrix() const { return priceMatrix_; }
        double getOptionPrice() const;
    };

} // namespace m2

#endif // CRANKNICOLSON_H
