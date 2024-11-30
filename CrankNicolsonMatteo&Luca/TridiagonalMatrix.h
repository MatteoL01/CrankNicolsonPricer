#ifndef TRIDIAGONALMATRIX_H
#define TRIDIAGONALMATRIX_H

#include <vector>
#include <stdexcept>

namespace m2 {

    class TridiagonalMatrix {
    private:
        std::vector<double> lower_;  // Subdiagonal (size n-1)
        std::vector<double> main_;   // Main diagonal (size n)
        std::vector<double> upper_;  // Superdiagonal (size n-1)
        size_t size_;                // Number of rows/columns

    public:
        explicit TridiagonalMatrix(size_t size, double initValue = 0.0)
            : lower_(size - 1, initValue), main_(size, initValue), upper_(size - 1, initValue), size_(size) {
        }

        size_t getSize() const { return size_; }

        double& lower(size_t i) {
            if (i >= size_ - 1) throw std::out_of_range("Lower diagonal index out of bounds");
            return lower_[i];
        }

        double& main(size_t i) {
            if (i >= size_) throw std::out_of_range("Main diagonal index out of bounds");
            return main_[i];
        }

        double& upper(size_t i) {
            if (i >= size_ - 1) throw std::out_of_range("Upper diagonal index out of bounds");
            return upper_[i];
        }

        const double& lower(size_t i) const { return lower_[i]; }
        const double& main(size_t i) const { return main_[i]; }
        const double& upper(size_t i) const { return upper_[i]; }

        // Thomas algorithm for solving Ax = b
        std::vector<double> solve(const std::vector<double>& rhs) const {
            if (rhs.size() != size_) throw std::invalid_argument("RHS vector size mismatch");

            std::vector<double> c_prime(size_ - 1);
            std::vector<double> d_prime(size_);
            std::vector<double> result(size_);

            // Forward sweep
            c_prime[0] = upper_[0] / main_[0];
            d_prime[0] = rhs[0] / main_[0];

            for (size_t i = 1; i < size_; ++i) {
                double denom = main_[i] - lower_[i - 1] * c_prime[i - 1];
                c_prime[i - 1] = upper_[i - 1] / denom;
                d_prime[i] = (rhs[i] - lower_[i - 1] * d_prime[i - 1]) / denom;
            }

            // Backward substitution
            result[size_ - 1] = d_prime[size_ - 1];
            for (int i = static_cast<int>(size_) - 2; i >= 0; --i) {
                result[i] = d_prime[i] - c_prime[i] * result[i + 1];
            }


            return result;
        }
    };

} // namespace m2

#endif // TRIDIAGONALMATRIX_H
