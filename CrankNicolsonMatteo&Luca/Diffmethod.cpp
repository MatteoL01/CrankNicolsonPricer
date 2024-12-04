#include "Diffmethod.h"

namespace m2
{
    American::American(Option& opt) : Option(opt), price_(0.0), T0prices_(M_ + 1, 0.0), delta_(M_, 0.0), gamma_(M_, 0.0), theta_(M_, 0.0), vega_(0.0), rho_(0.0), optionPrice_(0.0), boundary_(M_, 0.0)
    {
        dt_ = T_ / N_;
        opt.getCallPut() ? Smax_ = S0_ * 4 : Smax_ = K_ * 2;
        ds_ = Smax_ / M_;

        // initialize the matrix of price and greeks
        values_ = Matrix(N_ + 1, M_ + 1);

    }

    European::European(Option& opt) : Option(opt), price_(0.0), T0prices_(M_ + 1, 0.0), delta_(M_, 0.0), gamma_(M_, 0.0), theta_(M_, 0.0), vega_(0.0), rho_(0.0), optionPrice_(0.0), boundary_(M_, 0.0)
    {
        dt_ = T_ / N_;
        opt.getCallPut() ? Smax_ = S0_ * 4 : Smax_ = K_ * 2;
        ds_ = Smax_ / M_;

        // initialize the matrix of price and greeks
        values_ = Matrix(N_ + 1, M_ + 1);

    }


    void American::pricePut() {

        std::vector<double> a(M_ - 1), b(M_ - 1), c(M_ - 1), d(M_ - 1);

        for (unsigned int i = 0; i < (N_ + 1); i++)
        {
            // boundary condition asymptotic assumptions
            values_(i, 0) = K_; // when S = 0 the values is equal to K
            values_(i, M_) = 0; // when S = Smax the value is equal to zero
        }


        // create the two tridiagonal matrices
        Matrix T1(M_ - 1, M_ - 1), T2(M_ - 1, M_ - 1);

        // create vector of V and k (boundary times)
        std::vector<double> V(M_ - 1), k(M_ - 1);

        // Fill k vector
        k[0] = a[0] * Smax_;
        for (unsigned int i = 1; i < (M_ - 1); i++)
        {
            k[i] = 0.0;
        }

        // Add terminal values of V (at time N)
        for (unsigned int j = 0; j < (M_ - 1); j++)
        {
            V[j] = max(K_ - ds_ * (j + 1), 0);
            values_(0, j + 1) = V[j];
        }


        // vector W = T1 * V + k
        std::vector<double> W(M_ - 1);


        for (int n = N_; n > 0; n--)
        {

            // interpolate the interest rate
            double current_rate = interpolateRate(n * dt_, rates_);

            for (unsigned int j = 0; j < (M_ - 1); j++) {

                a[j] = 0.25 * (j + 1) * dt_ * (pow(sigma_, 2) * (j + 1) - current_rate);
                b[j] = (1 - 0.5 * pow(sigma_ * (j + 1), 2) * dt_);
                c[j] = 0.25 * (j + 1) * dt_ * (pow(sigma_, 2) * (j + 1) + current_rate);
                d[j] = (1 + (current_rate + 0.5 * pow(sigma_ * (j + 1), 2)) * dt_);
            }

            for (unsigned int i = 0; i < (M_ - 1); i++)
            {
                T1(i, i) = b[i];
                T2(i, i) = d[i];
            }

            for (unsigned int i = 0; i < (M_ - 2); i++)
            {
                T1(i, i + 1) = c[i];
                T1(i + 1, i) = a[i + 1];

                T2(i, i + 1) = -c[i];
                T2(i + 1, i) = -a[i + 1];

            }

            W = (T1 * V);
            W += k;


            crout(T2, W, V, M_);

            // Calculate exercise boundary
            double bound = 0.0;  // Default if no early exercise
            for (unsigned int i = 0; i < M_ - 1; i++) {
                double intrinsic_value = max(K_ - ds_ * (i + 1), 0.0);
                if (V[i] == intrinsic_value) {
                    bound = ds_ * (i + 1);  // Stock price where early exercise occurs
                    break;  // Only need the first occurrence
                }
            }
            boundary_[N_ - n] = bound;

            for (unsigned int i = 0; i < M_ - 1; i++)
            {
                V[i] = max(V[i], K_ - ds_ * (i + 1)); // american put payoff
                values_(N_ - n + 1, i + 1) = V[i];
            }

        }

        unsigned int pos = S0_ / ds_;
        price_ = values_(N_, pos);
        //printMatrix();

        T0prices_ = V;
    }

    void American::priceCall() {

        std::vector<double> a(M_ - 1), b(M_ - 1), c(M_ - 1), d(M_ - 1);

        for (unsigned int i = 0; i < (N_ + 1); i++)
        {
            // boundary condition asymptotic assumptions
            values_(i, 0) = 0; // when S = 0 the values is equal to 0
            values_(i, M_) = Smax_ - K_; // when S = Smax the value is equal to S - K
        }


        // create the two tridiagonal matrices
        Matrix T1(M_ - 1, M_ - 1), T2(M_ - 1, M_ - 1);

        // create vector of V and k (boundary times)
        std::vector<double> V(M_ - 1), k(M_ - 1);

        // Fill k vector
        for (unsigned int i = 0; i < (M_ - 2); i++)
        {
            k[i] = 0.0;
        }
        k[M_ - 2] = c[M_ - 2] * Smax_;

        // Add terminal values of V (at time N)
        for (unsigned int j = 0; j < (M_ - 1); j++)
        {
            V[j] = max(ds_ * (j + 1) - K_, 0);
            values_(0, j + 1) = V[j];
        }

        // vector W = T1 * V + k
        std::vector<double> W(M_ - 1);


        for (int n = N_; n > 0; n--)
        {

            // interpolate the interest rate
            double current_rate = interpolateRate(n * dt_, rates_);

            for (unsigned int j = 0; j < (M_ - 1); j++) {

                a[j] = 0.25 * (j + 1) * dt_ * (pow(sigma_, 2) * (j + 1) - current_rate);
                b[j] = (1 - 0.5 * pow(sigma_ * (j + 1), 2) * dt_);
                c[j] = 0.25 * (j + 1) * dt_ * (pow(sigma_, 2) * (j + 1) + current_rate);
                d[j] = (1 + (current_rate + 0.5 * pow(sigma_ * (j + 1), 2)) * dt_);
            }

            for (unsigned int i = 0; i < (M_ - 1); i++)
            {
                T1(i, i) = b[i];
                T2(i, i) = d[i];
            }

            for (unsigned int i = 0; i < (M_ - 2); i++)
            {
                T1(i, i + 1) = c[i];
                T1(i + 1, i) = a[i + 1];

                T2(i, i + 1) = -c[i];
                T2(i + 1, i) = -a[i + 1];

            }

            W = (T1 * V);
            W += k;


            crout(T2, W, V, M_);

            // Calculate exercise boundary
            double bound = 0.0;  // Default if no early exercise
            for (unsigned int i = 0; i < M_ - 1; i++) {
                double intrinsic_value = max( ds_ * (i + 1) - K_, 0.0);
                if (V[i] == intrinsic_value) {
                    bound = ds_ * (i + 1);  // Stock price where early exercise occurs
                    break;  // Only need the first occurrence
                }
            }
            boundary_[N_ - n] = bound;

            for (unsigned int i = 0; i < M_ - 1; i++)
            {
                V[i] = max(V[i], ds_ * (i + 1) - K_); // american call payoff 
                if (fabs(V[i]) < 1e-4) V[i] = 0.0;
                values_(N_ - n + 1, i + 1) = V[i];
            }

        }

        unsigned int pos = S0_ / ds_;
        price_ = values_(N_, pos);

        T0prices_ = V;
    }

    void European::pricePut()
    {

        std::vector<double> a(M_ - 1), b(M_ - 1), c(M_ - 1), d(M_ - 1);

        for (unsigned int i = 0; i < (N_ + 1); i++)
        {
            // boundary condition asymptotic assumptions
            values_(i, 0) = K_ * exp(-computeAverageRate(rates_, T_) * (T_ - dt_ * (N_ - i))); // when S = 0 the values is equal to K * e^(-rT) we consider r as the mean interest rate
            values_(i, M_) = 0; // when S = Smax the value is equal to zero
        }


        // create the two tridiagonal matrices
        Matrix T1(M_ - 1, M_ - 1), T2(M_ - 1, M_ - 1);

        // create vector of V and k (boundary times)
        std::vector<double> V(M_ - 1), k(M_ - 1);

        // Fill k vector
        k[0] = a[0] * Smax_; // check for this
        for (unsigned int i = 1; i < (M_ - 1); i++)
        {
            k[i] = 0.0;
        }

        // Add terminal values of V (at time N)
        for (unsigned int j = 0; j < (M_ - 1); j++)
        {
            V[j] = max(K_ - ds_ * (j + 1), 0);
            values_(0, j + 1) = V[j];
        }


        // vector W = T1 * V + k
        std::vector<double> W(M_ - 1);


        for (int n = N_; n > 0; n--)
        {

            // interpolate the interest rate
            double current_rate = interpolateRate(n * dt_, rates_);

            for (unsigned int j = 0; j < (M_ - 1); j++) {

                a[j] = 0.25 * (j + 1) * dt_ * (pow(sigma_, 2) * (j + 1) - current_rate);
                b[j] = (1 - 0.5 * pow(sigma_ * (j + 1), 2) * dt_);
                c[j] = 0.25 * (j + 1) * dt_ * (pow(sigma_, 2) * (j + 1) + current_rate);
                d[j] = (1 + (current_rate + 0.5 * pow(sigma_ * (j + 1), 2)) * dt_);
            }

            for (unsigned int i = 0; i < (M_ - 1); i++)
            {
                T1(i, i) = b[i];
                T2(i, i) = d[i];
            }

            for (unsigned int i = 0; i < (M_ - 2); i++)
            {
                T1(i, i + 1) = c[i];
                T1(i + 1, i) = a[i + 1];

                T2(i, i + 1) = -c[i];
                T2(i + 1, i) = -a[i + 1];

            }

            W = (T1 * V);
            W += k;


            crout(T2, W, V, M_);

            for (unsigned int i = 0; i < M_ - 1; i++)
            {
                V[i] = max(V[i], 0); // European payoff (consider only positive values of V[i] as it can get negative due to approximations
                values_(N_ - n + 1, i + 1) = V[i];
            }

        }

        unsigned int pos = S0_ / ds_;
        price_ = values_(N_, pos);
        //printMatrix();

        T0prices_ = V;
    }

    void European::priceCall()
    {
        std::vector<double> a(M_), b(M_), c(M_), d(M_);

        for (unsigned int i = 0; i < (N_ + 1); i++)
        {
            // boundary condition asymptotic assumptions
            values_(i, 0) = 0; // when S = 0 the values is equal to K
            values_(i, M_) = Smax_ - K_ * exp(-computeAverageRate(rates_, T_) * (T_ - dt_ * (N_ - i))); // when S = Smax the value is equal to [S - K * e^(-r(T-t))]
        }


        // create the two tridiagonal matrices
        Matrix T1(M_, M_), T2(M_, M_);

        // create vector of V and k (boundary times)
        std::vector<double> V(M_, M_), k(M_, M_);

        // Fill k vector
        for (unsigned int i = 0; i < M_; i++)
        {
            k[i] = 0.0;
        }
        k[M_] = c[M_] * Smax_; // check for this

        // Add terminal values of V (at time N)
        for (unsigned int j = 0; j < M_; j++)
        {
            V[j] = max(ds_ * (j + 1) - K_, 0);
            values_(0, j + 1) = V[j];
        }

        // vector W = T1 * V + k
        std::vector<double> W(M_);


        for (int n = N_; n > 0; n--)
        {
            k[M_] = c[M_] * (Smax_ - K_ * exp(-computeAverageRate(rates_, T_) * (T_ - dt_ * (N_ - n)))); // check for this
            // interpolate the interest rate
            double current_rate = interpolateRate(n * dt_, rates_);

            for (unsigned int j = 0; j < M_; j++) {

                a[j] = 0.25 * (j + 1) * dt_ * (pow(sigma_, 2) * (j + 1) - current_rate);
                b[j] = (1 - 0.5 * pow(sigma_ * (j + 1), 2) * dt_);
                c[j] = 0.25 * (j + 1) * dt_ * (pow(sigma_, 2) * (j + 1) + current_rate);
                d[j] = (1 + (current_rate + 0.5 * pow(sigma_ * (j + 1), 2)) * dt_);
            }

            for (unsigned int i = 0; i < M_; i++)
            {
                T1(i, i) = b[i];
                T2(i, i) = d[i];
            }

            for (unsigned int i = 0; i < M_ - 1; i++)
            {
                T1(i, i + 1) = c[i];
                T1(i + 1, i) = a[i + 1];

                T2(i, i + 1) = -c[i];
                T2(i + 1, i) = -a[i + 1];

            }

            W = (T1 * V);
            W += k;

            crout(T2, W, V, M_);

            for (unsigned int i = 0; i < M_; i++)
            {
                //V[i] = max(V[i],0); // european call payoff
                if (fabs(V[i]) < 1e-4) V[i] = 0.0;
                values_(N_ - n + 1, i + 1) = max(V[i], 0);
            }

        }

        unsigned int pos = S0_ / ds_;
        price_ = values_(N_, pos);
        //printMatrix();

        T0prices_ = V;

        
    }

    void European::calculateDelta()
    {
        for (unsigned int i = 0; i < M_; i++)
        {
            delta_[i] = (values_(N_, i + 1) - values_(N_, i)) / ds_;
        }
    }

    void American::calculateDelta()
    {
        for (unsigned int i = 0; i < M_; i++)
        {
            delta_[i] = (values_(N_, i + 1) - values_(N_, i)) / ds_;
        }
    }

    void European::calculateGamma()
    {
        for (unsigned int i = 0; i < M_; i++)
        {
            gamma_[i] = (delta_[i + 1] - delta_[i]) / ds_;
        }
    }

    void American::calculateGamma()
    {
        for (unsigned int i = 0; i < M_; i++)
        {
            gamma_[i] = (delta_[i + 1] - delta_[i]) / ds_;
        }
    }

    void European::calculateTheta()
    {
        for (unsigned int i = 0; i < M_; i++)
        {
            theta_[i] = (values_(N_ - 1, i) - values_(N_, i)) / dt_;
        }
    }

    void American::calculateTheta()
    {
        for (unsigned int i = 0; i < M_; i++)
        {
            theta_[i] = (values_(N_ -1, i) - values_(N_, i)) / dt_;
        }
    }

    void American::calculateVega()
    {
        optionPrice_ = price_;
        double epsilon = 0.001;

        sigma_ += epsilon;

        if (call_) { priceCall(); }
        else { pricePut(); }

        vega_ = (values_(N_, S0_ / ds_) - optionPrice_)/epsilon;

        sigma_ -= epsilon;
    }

    void European::calculateVega()
    {
        optionPrice_ = price_;
        double epsilon = 0.001;

        sigma_ += epsilon;
        
        if (call_) { priceCall(); }
        else { pricePut(); }

        vega_ = (values_(N_, S0_ / ds_) - optionPrice_) / epsilon;

        sigma_ -= epsilon;
    }

    void European::calculateRho()
    {
        
        double epsilon = 0.001;

        for (size_t i = 0; i < rates_.size() - 1; i++)
        {
            rates_[i].second += epsilon;
        }

        if (call_) { priceCall(); }
        else { pricePut(); }

        //std::cout << "Price: " << values_(N_, S0_ / ds_) - optionPrice_ << std::endl;

        rho_ = (values_(N_, S0_ / ds_) - optionPrice_) / epsilon;
    }

    void American::calculateRho()
    {
        double epsilon = 0.001;

        for (size_t i = 0; i < rates_.size() - 1; i++)
        {
            rates_[i].second += epsilon;
        }

        if (call_) { priceCall(); }
        else { pricePut(); }

        rho_ = (values_(N_, S0_ / ds_) - optionPrice_) / epsilon;
    }

}