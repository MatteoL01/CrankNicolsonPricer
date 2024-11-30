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

        for (unsigned int j = 0; j < (M_); j++) {
            float current_rate = interpolateRate(j * dt, rates_);

            a[j] = 0.25 * (j + 1) * dt * (pow(sigma_, 2) * (j + 1) - current_rate);
            b[j] = (1 - 0.5 * pow(ds * (j + 1), 2) * dt);
            c[j] = 0.25 * (j + 1) * dt * (pow(ds, 2) * (j + 1) + current_rate);
            d[j] = (1 + (current_rate + 0.5 * pow(ds * (j + 1), 2)) * dt);
        }

        // create the two tridiagonal matrices
        Matrix T1(M_, M_), T2(M_, M_);

        // create vector of V and k (boundary times)
        std::vector<float> V(M_), k(M_);

        // Add terminal values of V (at time N)
        for (unsigned int j = 0; j < M_; j++)
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