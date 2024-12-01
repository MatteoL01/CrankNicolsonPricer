#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>

namespace m2 {

    class Matrix
    {
    protected:
        int nl_;  // number of rows
        int nc_;  // number of columns
        std::vector<double> data_;  // stores matrix elements in a single 1D array

    public:
        // Constructors
        Matrix();
        Matrix(unsigned int nl, unsigned int nc);

        // Destructor
        ~Matrix();

        // Copy constructor
        Matrix(const Matrix& m1);

        // Getter methods
        int get_nl() const;
        int get_nc() const;

        // Access operators
        double operator()(unsigned int i, unsigned int j) const;
        double& operator()(unsigned int i, unsigned int j);

        // Assignment operator
        Matrix& operator=(const Matrix& m1);

    };

    // Stream output operator
    std::ostream& operator<<(std::ostream& st, const Matrix& m1);

    // Matrix arithmetic operators
    Matrix operator+(const Matrix& m1, const Matrix& m2);
    Matrix operator*(const Matrix& m1, const Matrix& m2);
    Matrix operator*(const Matrix& m1, const double d);
    Matrix operator*(const double d, const Matrix& m1);

    // Matrix-vector multiplication operator
    std::vector<double> operator*(const Matrix& mat, const std::vector<double>& vec);

} // namespace m2

#endif // MATRIX_H
