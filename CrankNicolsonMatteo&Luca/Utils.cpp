#include "Utils.h"

namespace m2 {
    double normalCDF(double x) {
        return 0.5 * erfc(-x * M_SQRT1_2);
    }

    double blackScholesPrice(bool isCall, double S, double K, double T, double r, double sigma) {
        double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
        double d2 = d1 - sigma * std::sqrt(T);

        if (isCall) {
            std::cout << "B-S price: " << S * normalCDF(d1) - K * std::exp(-r * T) * normalCDF(d2) << std::endl;
            std::cout << "Delta: " << normalCDF(d1) << std::endl;
            return S * normalCDF(d1) - K * std::exp(-r * T) * normalCDF(d2);
        }
        else {
            std::cout << "B-S price: " << K * std::exp(-r * T) * normalCDF(-d2) - S * normalCDF(-d1) << std::endl;
            std::cout << "Delta: " << normalCDF(d1) - 1 << std::endl;
            return K * std::exp(-r * T) * normalCDF(-d2) - S * normalCDF(-d1);
        }

        
    }

    float max(float a, float b) {
        return (a > b) ? a : b;
    }

    float interpolateRate(double t, const std::vector<std::pair<double, double>>& rates) {
        if (t < rates.front().first) {
            return rates.front().second; // Clamp to the first rate
        }
        if (t > rates.back().first) {
            return rates.back().second; // Clamp to the last rate
        }
        for (size_t i = 0; i < rates.size() - 1; i++) {
            if (t >= rates[i].first && t <= rates[i + 1].first) {
                float t1 = rates[i].first;
                float t2 = rates[i + 1].first;
                float r1 = rates[i].second;
                float r2 = rates[i + 1].second;
                return r1 + (t - t1) / (t2 - t1) * (r2 - r1); // Linear interpolation
            }
        }
        return rates.back().second; // Redundant but ensures completeness
    }

    std::vector<double>& operator+=(std::vector<double>& lhs, const std::vector<double>& rhs) {
        if (lhs.size() != rhs.size()) {
            throw std::invalid_argument("Vectors must have the same size for addition.");
        }

        for (size_t i = 0; i < lhs.size(); ++i) {
            lhs[i] += rhs[i];
        }
        return lhs;
    }

    void crout(Matrix& T2, std::vector<double>& W, std::vector<double>& V, int M) {
        M = M - 1; // we want dimension M - 1

        Matrix L(M, M), U(M, M);

        std::vector<float> Z(M);

        L(0, 0) = T2(0, 0);
        U(0, 1) = T2(0, 0) / L(0, 0);
        for (unsigned int i = 1; i <= (M - 2); i++) {
            L(i, i - 1) = T2(i, i - 1);
            L(i, i) = T2(i, i) - L(i, i - 1) * U(i - 1, i);
            U(i, i + 1) = T2(i, i + 1) / L(i, i);
        }

        L(M - 1, M - 2) = T2(M - 1, M - 2);
        L(M - 1, M - 1) = T2(M - 1, M - 1) - L(M - 1, M - 2) * U(M - 2, M - 1);

        // solve the system
        // first solve Lz = b for z, where z = U V
        Z[0] = W[0];
        for (unsigned int i = 1; i <= (M - 1); i++)
        {
            Z[i] = (W[i] - L(i, i - 1) * Z[i - 1]) / L(i, i);
        }

        // solve U V = Z
        V[M - 1] = Z[M - 1];
        for (int i = (M - 2); i >= 0; i--)
        {
            V[i] = Z[i] - U(i, i + 1) * V[i + 1];
        }


    }

    double computeAverageRate(const std::vector<std::pair<double, double>>& rates, double T) {
        if (rates.empty())
            throw std::runtime_error("Rates vector is empty");

        if (T <= rates.front().first)
            return rates.front().second; // Return the first rate for times before the first timestamp

        double totalRate = 0.0;
        double totalWeight = 0.0;

        for (size_t i = 1; i < rates.size(); ++i) {
            double start = rates[i - 1].first;
            double end = rates[i].first;
            double rate = rates[i - 1].second;

            if (T <= start)
                break; // No contribution from intervals completely after T

            double timeDiff = std::min(end, T) - start; // Consider up to T if within this interval
            totalRate += rate * timeDiff;
            totalWeight += timeDiff;

            if (T <= end)
                break; // Stop after reaching T
        }

        // Handle the remaining time beyond the last known rate
        if (T > rates.back().first) {
            double timeDiff = T - rates.back().first;
            totalRate += rates.back().second * timeDiff;
            totalWeight += timeDiff;
        }

        if (totalWeight == 0.0)
            throw std::runtime_error("Total weight is zero, invalid time intervals");

        return totalRate / totalWeight;
    }

}