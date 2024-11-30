#ifndef PRICER_H
#define PRICER_H

#include <vector>
#include <utility>
#include <string>
#include "Option.h"

namespace m2 {
    class VanillaOption : public Option  // Inherit from Option class
    {
    protected:
        Discretization discretization_;

    public:
        VanillaOption(const Discretization& discretization, const std::vector<std::pair<double, double>>& rates)
            : Option(), discretization_(discretization)
        {
        }

        virtual ~VanillaOption() = default;

        double get_interest_rate(double time) const;

        int get_M() const { return discretization_.get_M(); }
        int get_N() const { return discretization_.get_N(); }
        double get_h() const { return discretization_.get_h(); }
        double get_dt() const { return discretization_.get_dt(); }

        virtual void do_pricing_and_create_files_text() = 0;
    };

    class EuropeanOption : public VanillaOption
    {
    protected:
        int put_call_; // 1 for a put, -1 for a call
        int ind_s_0_;  // index i such that S[i] >= S_0

    public:
        EuropeanOption(int N, double x_max, int M, double T, double K, double S_0, double sigma, const std::vector<std::pair<double, double>>& rates, int put_call)
            : VanillaOption(Discretization(N, x_max, M, T), rates), put_call_(put_call) {
        }

        void do_pricing_and_create_files_text();
        void fill_P(std::vector<float>& A_a, std::vector<float>& A_b, std::vector<float>& A_c,
            std::vector<float>& B_a, std::vector<float>& B_b, std::vector<float>& B_c,
            const std::vector<float>& S, std::vector<float>& P);
        void compute(std::vector<float>& A_a, std::vector<float>& A_b, std::vector<float>& A_c,
            std::vector<float>& B_a, std::vector<float>& B_b, std::vector<float>& B_c,
            std::vector<float>& time, std::vector<float>& P, float& prev_price);
        void create_files_txt(std::vector<float>& A_a, std::vector<float>& A_b, std::vector<float>& A_c,
            std::vector<float>& B_a, std::vector<float>& B_b, std::vector<float>& B_c,
            std::vector<float>& S, std::vector<float>& time, std::vector<float>& P, float& prev_price);
        void create_file_txt_greeks(std::vector<float>& A_a, std::vector<float>& A_b, std::vector<float>& A_c,
            std::vector<float>& B_a, std::vector<float>& B_b, std::vector<float>& B_c,
            std::vector<float>& S, std::vector<float>& time, std::vector<float>& P, float& prev_price);
    };

    class AmericanOption : public VanillaOption
    {
    protected:
        int put_call_; // 1 for a put, -1 for a call
        int ind_s_0_;  // index i such that S[i] >= S_0

    public:
        AmericanOption(int N, double x_max, int M, double T, double K, double S_0, double sigma, const std::vector<std::pair<double, double>>& rates, int put_call)
            : VanillaOption(Discretization(N, x_max, M, T), rates), put_call_(put_call) {
        }

        void do_pricing_and_create_files_text();
        void fill_P(std::vector<float>& A_a, std::vector<float>& A_b, std::vector<float>& A_c,
            std::vector<float>& B_a, std::vector<float>& B_b, std::vector<float>& B_c,
            std::vector<float>& S, std::vector<float>& P);
        void compute(std::vector<float>& A_a, std::vector<float>& A_b, std::vector<float>& A_c,
            std::vector<float>& B_a, std::vector<float>& B_b, std::vector<float>& B_c,
            std::vector<float>& S, std::vector<float>& time, std::vector<float>& P,
            std::vector<float>& Frontier, float& prev_price);
        void create_files_txt(std::vector<float>& A_a, std::vector<float>& A_b, std::vector<float>& A_c,
            std::vector<float>& B_a, std::vector<float>& B_b, std::vector<float>& B_c,
            std::vector<float>& S, std::vector<float>& time, std::vector<float>& P,
            std::vector<float>& Frontier, float& prev_price);
        void create_file_txt_greeks(std::vector<float>& A_a, std::vector<float>& A_b, std::vector<float>& A_c,
            std::vector<float>& B_a, std::vector<float>& B_b, std::vector<float>& B_c,
            std::vector<float>& S, std::vector<float>& time, std::vector<float>& P,
            std::vector<float>& Frontier, float& prev_price);
    };
}

#endif // PRICER_H
