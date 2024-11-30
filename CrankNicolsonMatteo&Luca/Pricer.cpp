#include "pricer.h"
#include "Utils.h"
#include <fstream>
#include <iostream>
#include <cmath>

// VanillaOption
m2::VanillaOption::VanillaOption(const Discretization& discretization, const std::vector<std::pair<double, double>>& rates)
    : Option(), discretization_(discretization) {
    // Initialization here (if necessary)
}

double m2::VanillaOption::get_interest_rate(double time) const {
    // Ensure that the time is within the bounds of the rates vector
    if (time <= rates_.front().first) {
        return rates_.front().second; // Return the first rate if time is before the first interval
    }

    for (size_t i = 0; i < rates_.size() - 1; ++i) {
        if (time >= rates_[i].first && time < rates_[i + 1].first) {
            // Perform linear interpolation between the two points
            double t0 = rates_[i].first;
            double r0 = rates_[i].second;
            double t1 = rates_[i + 1].first;
            double r1 = rates_[i + 1].second;

            // Linear interpolation formula
            double rate = r0 + (time - t0) * (r1 - r0) / (t1 - t0);
            return rate;
        }
    }

    // If time exceeds the last interval, return the last rate
    return rates_.back().second;
}


// EuropeanOption
m2::EuropeanOption::EuropeanOption(int N, double x_max, int M, double T, double K, double S_0, double sigma, const std::vector<std::pair<double, double>>& rates, int put_call)
    : VanillaOption(Discretization(N, x_max, M, T), rates), put_call_(put_call) {
    // Initialization here (if needed)
}

void m2::EuropeanOption::do_pricing_and_create_files_text() {
    // Create time discretization set
    std::vector<float> time(get_M() + 1);
    for (int m = 0; m <= get_M(); ++m) {
        time[m] = m * get_dt();
    }

    // Create spot discretization set
    std::vector<float> S(get_N() + 1);
    for (int j = 0; j <= get_N(); ++j) {
        S[j] = j * get_h();
    }

    // Set field ind_s_0_
    ind_s_0_ = 0;
    while (S[ind_s_0_] < getS0()) {
        ++ind_s_0_;
        if (ind_s_0_ >= S.size()) break;  // Check if it goes out of bounds
    }
    ind_s_0_ = std::max(ind_s_0_ - 1, 0);  // Prevent going out of bounds

    std::vector<float> A_a(get_N());
    std::vector<float> A_b(get_N() + 1);
    std::vector<float> A_c(get_N());
    std::vector<float> B_a(get_N());
    std::vector<float> B_b(get_N() + 1);
    std::vector<float> B_c(get_N());
    std::vector<float> P(get_N() + 1);
    float prev_price = 0.0f;

    // Fill the price vector
    fill_P(A_a, A_b, A_c, B_a, B_b, B_c, S, P);

    // Compute the solution prices
    compute(A_a, A_b, A_c, B_a, B_b, B_c, time, P, prev_price);

    // Create .txt files for outputs of P, delta, and other Greeks
    create_files_txt(A_a, A_b, A_c, B_a, B_b, B_c, S, time, P, prev_price);
}


void m2::EuropeanOption::fill_P(std::vector<float>& A_a, std::vector<float>& A_b, std::vector<float>& A_c,
    std::vector<float>& B_a, std::vector<float>& B_b, std::vector<float>& B_c,
    const std::vector<float>& S, std::vector<float>& P)
{
    // Initialize the A and B coefficient vectors
    init_A(A_a, A_b, A_c, S, get_N(), rates_, sigma_, get_dt(), get_h());
    init_B(B_a, B_b, B_c, S, get_N(), rates_, sigma_, get_dt(), get_h());

    // Fill the price vector P with the option's payoff (European put/call)
    for (int j = 0; j < get_N() + 1; ++j) {
        P[j] = std::max(float(put_call_ * (K_ - S[j])), 0.0f);  // Calculate the payoff (put or call option)
    }
}


void m2::EuropeanOption::compute(std::vector<float>& A_a, std::vector<float>& A_b, std::vector<float>& A_c,
    std::vector<float>& B_a, std::vector<float>& B_b, std::vector<float>& B_c,
    std::vector<float>& time, std::vector<float>& P, float& prev_price)
{
    for (int m = get_M(); m > -1; m--) {
        // Calculate the put/call option payoff at the boundary
        P[int(float(get_N()) / 2 * (1 - put_call_))] = K_ * exp(-r_ * put_call_ * (T_ - time[m]));
        P[int(float(get_N()) / 2 * (put_call_ + 1))] = 0.;

        // Store the option price for m = 0
        if (m == 0) { prev_price = P[ind_s_0_]; }

        // Perform matrix multiplication for the tridiagonal system
        matmul_tridiag(P, B_a, B_b, B_c, get_N() + 1);

        // Solve the tridiagonal system using the Thomas algorithm
        thomas_alg(P, A_a, A_b, A_c, P, get_N() + 1);
    }
}


void m2::EuropeanOption::create_files_txt(std::vector<float>& A_a, std::vector<float>& A_b, std::vector<float>& A_c,
    std::vector<float>& B_a, std::vector<float>& B_b, std::vector<float>& B_c,
    std::vector<float>& S, std::vector<float>& time, std::vector<float>& P, float& prev_price)
{
    // Output for P as a function of S
    std::fstream file_P("P_output", std::ios::out);
    for (int j = 0; j < get_N() + 1; j++) {
        file_P << P[j] << "\n";
    }

    // Output for delta as a function of S
    std::fstream file_D("Delta_output", std::ios::out);
    for (int j = 1; j < get_N(); j++) {
        file_D << (P[j + 1] - P[j - 1]) / (2 * get_h()) << "\n";
    }

    // Output for the greeks
    create_file_txt_greeks(A_a, A_b, A_c, B_a, B_b, B_c, S, time, P, prev_price);
}


// AmericanOption
m2::AmericanOption::AmericanOption(int N, double x_max, int M, double T, double K, double S_0, double sigma, const std::vector<std::pair<double, double>>& rates, int put_call)
    : VanillaOption(Discretization(N, x_max, M, T), rates), put_call_(put_call) {
    // Initialization here (if needed)
}

void m2::AmericanOption::do_pricing_and_create_files_text() {
    // Implement pricing and text file creation for American options
    // Placeholder logic:
    float* A_a = new float[get_M()];
    float* A_b = new float[get_M()];
    float* A_c = new float[get_M()];
    float* B_a = new float[get_M()];
    float* B_b = new float[get_M()];
    float* B_c = new float[get_M()];
    float* S = new float[get_M()];
    float* P = new float[get_M()];
    float* Frontier = new float[get_M()];
    float* time = new float[get_M()];

    float prev_price = 0.0;

    fill_P(A_a, A_b, A_c, B_a, B_b, B_c, S, P);
    compute(A_a, A_b, A_c, B_a, B_b, B_c, S, time, P, Frontier, prev_price);
    create_files_txt(A_a, A_b, A_c, B_a, B_b, B_c, S, time, P, Frontier, prev_price);

    // Clean up dynamically allocated memory
    delete[] A_a;
    delete[] A_b;
    delete[] A_c;
    delete[] B_a;
    delete[] B_b;
    delete[] B_c;
    delete[] S;
    delete[] P;
    delete[] Frontier;
    delete[] time;
}

void m2::AmericanOption::fill_P(float* A_a, float* A_b, float* A_c, float* B_a, float* B_b, float* B_c, float* S, float* P) {
    // Fill the matrices A and B, and the vector P for American option
    // Implement the logic for filling the matrices here
}

void m2::AmericanOption::compute(float* A_a, float* A_b, float* A_c, float* B_a, float* B_b, float* B_c, float* S, float* time, float* P, float* Frontier, float& prev_price) {
    // Compute the option price for the American option
    // This function will update P and prev_price and calculate the American option's early exercise frontier
}

void m2::AmericanOption::create_files_txt(float* A_a, float* A_b, float* A_c, float* B_a, float* B_b, float* B_c, float* S, float* time, float* P, float* Frontier, float& prev_price) {
    // Create a text file with the results of the American option pricing
    std::ofstream file("american_option.txt");
    if (file.is_open()) {
        file << "Option Price: " << prev_price << std::endl;
        // Output other data to the file as needed
        file.close();
    }
    else {
        std::cerr << "Unable to open file for writing!" << std::endl;
    }
}

