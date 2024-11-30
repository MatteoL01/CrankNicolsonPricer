#include "Diffmethod.h"

namespace m2
{
	American::American()
	{
		values_(timeDiscr_, spotDiscr_);
	}
	
	European::European()
	{
		values_(timeDiscr_, spotDiscr_);
	}

    void American::pricePut() {
        float Smax = 2 * K_;
        float dt = T_ / N_;
        float ds = Smax / M_;

        std::vector<float> a(M_), b(M_), c(M_), d(M_);

        for (unsigned int i = 0; i < (N_); i++)
        {   
            // boundary condition asymptotic assumptions
            values_(i, 0) = K_; // when S = 0 the values is equal to K
            values_(i, M_ - 1) = 0; // when S = Smax the value is equal to zero
        }

        // interpolate the interest rate
        float current_rate = interpolateRate(N_ * dt, rates_);

        for (unsigned int j = 0; j < (M_ - 1); j++) {

            a[j] = 0.25 * (j + 1) * dt * (pow(sigma_, 2) * (j + 1) - current_rate);
            b[j] = (1 - 0.5 * pow(ds * (j + 1), 2) * dt);
            c[j] = 0.25 * (j + 1) * dt * (pow(ds, 2) * (j + 1) + current_rate);
            d[j] = (1 + (current_rate + 0.5 * pow(ds * (j + 1), 2)) * dt);
        }

        // create the two tridiagonal matrices
        Matrix T1(M_ - 1, M_ - 1), T2(M_ - 1, M_ - 1);

        for (unsigned int i = 0; i < (M_ - 1); i++)
        {
            T1(i, i) = b[i];
            T2(i, i) = d[i];
        }

        for (unsigned int i = 0; i < (M_ - 2); i++)
        {
            T1(i, i +1) = c[i];
            T1(i + 1, i) = a[i + 1];

            T2(i, i + 1) = -c[i];
            T2(i + 1, i) = -a[i+1];

        }



        // create vector of V and k (boundary times)
        std::vector<float> V(M_), k(M_);

        // Add terminal values of V (at time N)
        for (unsigned int j = 0; j < M_-1; j++)
        {
            V[j] = max(K_ - ds * (j + 1), 0); 
            values_(0,j) = V[j];
        }

        // vector W = T1 * V + k
        std::vector<float> W(M_);

        W = (T1 * V);
        W += k;




    }

}