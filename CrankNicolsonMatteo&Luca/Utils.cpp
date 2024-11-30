#include "Utils.h"

double normalCDF(double x) {
    return 0.5 * erfc(-x * M_SQRT1_2);
}

double blackScholesPrice(bool isCall, double S, double K, double T, double r, double sigma) {
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);

    if (isCall) {
        return S * normalCDF(d1) - K * std::exp(-r * T) * normalCDF(d2);
    }
    else {
        return K * std::exp(-r * T) * normalCDF(-d2) - S * normalCDF(-d1);
    }
}

float max(float a, float b) {
    return (a > b) ? a : b;
}

float interpolateRate(float t, const std::vector<std::pair<float, float>>& rates) {
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

std::vector<float>& operator+=(std::vector<float>& lhs, const std::vector<float>& rhs) {
    if (lhs.size() != rhs.size()) {
        throw std::invalid_argument("Vectors must have the same size for addition.");
    }

    for (size_t i = 0; i < lhs.size(); ++i) {
        lhs[i] += rhs[i];
    }
    return lhs;
}

void crout(m2::Matrix& T2, std::vector<float>& W, std::vector<float>& V, m2::Matrix& result, unsigned int M) {
    m2::Matrix L(M, M), U(M, M);

    L(0, 0) = T2(0, 0);
    U(0, 1) = T2(0, 0) / L(0, 0);
    for (unsigned int i = 1; i <= (M - 2); i++) {
        L(i, i - 1) = T2(i, i - 1);
        L(i, i) = T2(i, i) - L(i, i - 1) * U(i - 1, i);
        U(i, i + 1) = T2(i, i + 1) / L(i, i);
    }

    L(M - 1, M - 2) = T2(M - 1, M - 2);
    L(M - 1, M - 1) = T2(M - 1, M - 1) - L(M - 1, M - 2) * U(M - 2, M - 1);
}

double computeAverageRate(const std::vector<std::pair<double, double>>& rates, double T) {
    if (rates.empty()) throw std::runtime_error("Rates vector is empty");

    double totalRate = 0.0;
    double totalWeight = 0.0;

    for (size_t i = 1; i < rates.size(); ++i) {
        double timeDiff = rates[i].first - rates[i - 1].first;
        totalRate += rates[i - 1].second * timeDiff;
        totalWeight += timeDiff;
    }

    // Account for time remaining till T
    double lastTime = rates.back().first;
    if (T > lastTime) {
        totalRate += rates.back().second * (T - lastTime);
        totalWeight += (T - lastTime);
    }

    return totalRate / totalWeight;
}
