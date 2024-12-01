#include "Matrix.h"

namespace m2 {

    Matrix::Matrix(): nl_(1), nc_(1), data_(nl_* nc_, 0.0) {}

    // Constructor with dimensions, initializes with zero
    Matrix::Matrix(unsigned int nl, unsigned int nc)
        : nl_(nl), nc_(nc), data_(nl* nc, 0.0) {
    }


    // Destructor (handled by std::vector)
    Matrix::~Matrix() {}

    // Copy constructor (deep copy)
    Matrix::Matrix(const Matrix& m1)
        : nl_(m1.nl_), nc_(m1.nc_), data_(m1.data_) {
    }

    // Get number of rows
    int Matrix::get_nl() const { return nl_; }

    // Get number of columns
    int Matrix::get_nc() const { return nc_; }

    // Access operator (const version)
    float Matrix::operator()(unsigned int i, unsigned int j) const {
        if (i >= nl_ || j >= nc_ || i < 0 || j < 0) {
            throw std::out_of_range("Index out of bounds");
        }
        return data_[i * nc_ + j];
    }

    // Access operator (non-const version)
    float& Matrix::operator()(unsigned int i, unsigned int j) {
        if (i >= nl_ || j >= nc_ || i < 0 || j < 0) {
            throw std::out_of_range("Index out of bounds");
        }
        return data_[i * nc_ + j];
    }

    // Assignment operator (deep copy)
    Matrix& Matrix::operator=(const Matrix& m1) {
        if (this != &m1) {
            nl_ = m1.nl_;
            nc_ = m1.nc_;
            data_ = m1.data_;
        }
        return *this;
    }

    // Matrix extraction (this function could vary depending on the specific requirement)
    Matrix Matrix::extract(int index) {
        if (index < 0 || index >= nl_) {
            throw std::out_of_range("Invalid index");
        }
        Matrix result(1, nc_);
        for (int j = 0; j < nc_; ++j) {
            result(0, j) = (*this)(index, j);
        }
        return result;
    }

    // Output stream operator
    std::ostream& operator<<(std::ostream& st, const Matrix& m1) {
        for (int i = 0; i < m1.get_nl(); ++i) {
            for (int j = 0; j < m1.get_nc(); ++j) {
                st << std::setw(2) << std::setprecision(3) << m1(i, j) << " ";
            }
            st << std::endl;
        }
        return st;
    }

    std::vector<float> operator*(const Matrix& mat, const std::vector<float>& vec) {
        if (mat.get_nc() != vec.size()) {
            throw std::invalid_argument("Matrix columns must match vector size.");
        }

        std::vector<float> result(mat.get_nl(), 0.0f); // Initialize result vector with 0s

        for (int i = 0; i < mat.get_nl(); ++i) {
            for (int j = 0; j < mat.get_nc(); ++j) {
                result[i] += mat(i, j) * vec[j];
            }
        }

        return result;
    }

    // Matrix addition
    Matrix operator+(const Matrix& m1, const Matrix& m2) {
        if (m1.get_nl() != m2.get_nl() || m1.get_nc() != m2.get_nc()) {
            throw std::invalid_argument("Matrix dimensions must match for addition");
        }
        Matrix result(m1.get_nl(), m1.get_nc());
        for (int i = 0; i < m1.get_nl(); ++i) {
            for (int j = 0; j < m1.get_nc(); ++j) {
                result(i, j) = m1(i, j) + m2(i, j);
            }
        }
        return result;
    }

    // Matrix multiplication
    Matrix operator*(const Matrix& m1, const Matrix& m2) {
        if (m1.get_nc() != m2.get_nl()) {
            throw std::invalid_argument("Matrix dimensions must match for multiplication");
        }
        Matrix result(m1.get_nl(), m2.get_nc());
        for (int i = 0; i < m1.get_nl(); ++i) {
            for (int j = 0; j < m2.get_nc(); ++j) {
                double sum = 0;
                for (int k = 0; k < m1.get_nc(); ++k) {
                    sum += m1(i, k) * m2(k, j);
                }
                result(i, j) = sum;
            }
        }
        return result;
    }

    // Scalar multiplication
    Matrix operator*(const Matrix& m1, const double d) {
        Matrix result(m1.get_nl(), m1.get_nc());
        for (int i = 0; i < m1.get_nl(); ++i) {
            for (int j = 0; j < m1.get_nc(); ++j) {
                result(i, j) = m1(i, j) * d;
            }
        }
        return result;
    }

    // Scalar multiplication (commutative)
    Matrix operator*(const double d, const Matrix& m1) {
        return m1 * d;
    }

} // namespace m2
