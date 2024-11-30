#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <vector>
#include <utility> 
#include <algorithm>
#include <stdexcept>

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440 // Value of 1 / sqrt(2)
#endif

// Normal cumulative distribution function
double normalCDF(double x) {
    return 0.5 * erfc(-x * M_SQRT1_2);
}

// Black-Scholes formula for European call and put options
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

double getInterestRate(const std::vector<std::pair<double, double>>& rates, double time) {
    // Find the surrounding rates for interpolation
    for (size_t i = 0; i < rates.size() - 1; i++) {
        if (time >= rates[i].first && time <= rates[i + 1].first) {
            double t0 = rates[i].first;
            double r0 = rates[i].second;
            double t1 = rates[i + 1].first;
            double r1 = rates[i + 1].second;

            // Linear interpolation
            return r0 + (time - t0) * (r1 - r0) / (t1 - t0);
        }
    }
    // If time is outside the range, return the last rate (or handle this case differently)
    return rates.back().second;  // Use the last rate if out of range
}

// Compute the average rate over time using the piecewise rate structure
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

void thomas_alg(std::vector<float>& x, const std::vector<float>& a, const std::vector<float>& b, const std::vector<float>& c, const std::vector<float>& d, int n) {
    // Initialize w and g as vectors
    std::vector<float> w(n - 1);
    std::vector<float> g(n);

    // Forward sweep
    w[0] = c[0] / b[0];
    g[0] = d[0] / b[0];

    for (int i = 1; i < n - 1; i++) {
        w[i] = c[i] / (b[i] - a[i - 1] * w[i - 1]);
    }

    for (int i = 1; i < n; i++) {
        g[i] = (d[i] - a[i - 1] * g[i - 1]) / (b[i] - a[i - 1] * w[i - 1]);
    }

    // Backward substitution
    x[n - 1] = g[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = g[i] - w[i] * x[i + 1];
    }
}

void matmul_tridiag(std::vector<float>& x, const std::vector<float>& a, const std::vector<float>& b, const std::vector<float>& c, int n) {
    float x_2;
    float x_1 = x[0];
    x[0] = b[0] * x[0] + c[0] * x[1];

    for (int i = 1; i < n - 1; ++i) {
        x_2 = x[i];
        x[i] = a[i - 1] * x_1 + b[i] * x[i] + c[i] * x[i + 1];
        x_1 = x_2;
    }
    x[n - 1] = a[n - 2] * x_1 + b[n - 1] * x[n - 1];
}


void init_A(std::vector<float>& A_a, std::vector<float>& A_b, std::vector<float>& A_c, const std::vector<float>& S, int N, const std::vector<std::pair<double, double>>& rates, float sigma, float dt, float h) {
    for (int j = 1; j < N; ++j) {
        // Find the interest rate corresponding to the current time step j
        double r = getInterestRate(rates, j * dt);  // assuming time steps are indexed by dt

        A_a[j - 1] = (r * S[j]) / (4 * h) - (sigma * sigma * S[j] * S[j]) / (4 * h * h);
        A_b[j] = (1 / dt + r + (sigma * sigma * S[j] * S[j]) / (2 * h * h));
        A_c[j] = -(r * S[j]) / (4 * h) - (sigma * sigma * S[j] * S[j]) / (4 * h * h);
    }
    A_c[0] = 0.0f;
    A_a[N - 1] = 0.0f;
    A_b[0] = 1.0f;
    A_b[N] = 1.0f;
}


void init_B(std::vector<float>& B_a, std::vector<float>& B_b, std::vector<float>& B_c, const std::vector<float>& S, int N, const std::vector<std::pair<double, double>>& rates, float sigma, float dt, float h) {
    for (int j = 1; j < N; ++j) {
        // Find the interest rate corresponding to the current time step j
        double r = getInterestRate(rates, j * dt);  // assuming time steps are indexed by dt

        B_a[j - 1] = -(r * S[j]) / (4 * h) + (sigma * sigma * S[j] * S[j]) / (4 * h * h);
        B_b[j] = 1 / dt - (sigma * sigma * S[j] * S[j]) / (2 * h * h);
        B_c[j] = (r * S[j]) / (4 * h) + (sigma * sigma * S[j] * S[j]) / (4 * h * h);
    }
    B_c[0] = 0.0f;
    B_a[N - 1] = 0.0f;
    B_b[0] = 1.0f;
    B_b[N] = 1.0f;
}




#endif // UTILS_H
