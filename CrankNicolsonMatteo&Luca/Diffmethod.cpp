#include "Diffmethod.h"

namespace m2
{
	American::American(Option& opt) : call_(opt.getCallPut()), T_(opt.getT()), K_(opt.getK()), N_(opt.getTimeDiscr()), M_(opt.getSpotDiscr()),
        S0_(opt.getS0()), sigma_(opt.getSigma()), rates_(opt.getRates())
	{   
        values_ = Matrix(N_+1, M_+1);
	}
	
	European::European(Option& opt) : call_(opt.getCallPut()), T_(opt.getT()), K_(opt.getK()), N_(opt.getTimeDiscr()), M_(opt.getSpotDiscr()),
        S0_(opt.getS0()), sigma_(opt.getSigma()), rates_(opt.getRates())
	{
        values_ = Matrix(N_+1,M_+1);
	}

    void American::pricePut() {
        

        float Smax = 2 * K_;
        float dt = T_ / N_;
        float ds = Smax / M_;

        std::vector<float> a(M_ - 1), b(M_ - 1), c(M_ - 1), d(M_ - 1);

        for (unsigned int i = 0; i < (N_+1); i++)
        {   
            // boundary condition asymptotic assumptions
            values_(i, 0) = K_; // when S = 0 the values is equal to K
            values_(i, M_) = 0; // when S = Smax the value is equal to zero
        }

                
        // create the two tridiagonal matrices
        Matrix T1(M_ - 1, M_ - 1), T2(M_ - 1, M_ - 1);

        // create vector of V and k (boundary times)
        std::vector<float> V(M_-1), k(M_-1);

        // Fill k vector
        k[0] = a[0] * Smax;
        for (unsigned int i = 1; i < (M_ - 1); i++)
        {
            k[i] = 0.0;
        }
      
        // Add terminal values of V (at time N)
        for (unsigned int j = 0; j < (M_ - 1); j++)
        {
            V[j] = max(K_ - ds * (j + 1), 0); 
            values_(0,j+1) = V[j];
        }


        // vector W = T1 * V + k
        std::vector<float> W(M_-1);

        
        for (int n = N_; n > 0; n--)
        {
            
            // interpolate the interest rate
            float current_rate = interpolateRate(n * dt, rates_);

            for (unsigned int j = 0; j < (M_ - 1); j++) {

                a[j] = 0.25 * (j + 1) * dt * (pow(sigma_, 2) * (j + 1) - current_rate);
                b[j] = (1 - 0.5 * pow(sigma_ * (j + 1), 2) * dt);
                c[j] = 0.25 * (j + 1) * dt * (pow(sigma_, 2) * (j + 1) + current_rate);
                d[j] = (1 + (current_rate + 0.5 * pow(sigma_ * (j + 1), 2)) * dt);
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
                V[i] = max(V[i], K_ - ds * (i + 1)); // american put payoff
                values_(N_ - n + 1, i + 1) = V[i];
            }  

        }

        
        printMatrix();
    }

    void American::priceCall() {


        float Smax = 2 * S0_;
        float dt = T_ / N_;
        float ds = Smax / M_;

        std::vector<float> a(M_ - 1), b(M_ - 1), c(M_ - 1), d(M_ - 1);

        for (unsigned int i = 0; i < (N_ + 1); i++)
        {
            // boundary condition asymptotic assumptions
            values_(i, 0) = 0; // when S = 0 the values is equal to K
            values_(i, M_) = Smax - K_; // when S = Smax the value is equal to zero
        }


        // create the two tridiagonal matrices
        Matrix T1(M_ - 1, M_ - 1), T2(M_ - 1, M_ - 1);

        // create vector of V and k (boundary times)
        std::vector<float> V(M_ - 1), k(M_ - 1);

        // Fill k vector
        for (unsigned int i = 0; i < (M_ - 2); i++)
        {
            k[i] = 0.0;
        }
        k[M_ - 2] = c[M_ -2] * Smax;

        // Add terminal values of V (at time N)
        for (unsigned int j = 0; j < (M_ - 1); j++)
        {
            V[j] = max(ds * (j + 1) - K_, 0);
            values_(0, j + 1) = V[j];
        }

        // vector W = T1 * V + k
        std::vector<float> W(M_ - 1);


        for (int n = N_; n > 0; n--)
        {

            // interpolate the interest rate
            float current_rate = interpolateRate(n * dt, rates_);

            for (unsigned int j = 0; j < (M_ - 1); j++) {

                a[j] = 0.25 * (j + 1) * dt * (pow(sigma_, 2) * (j + 1) - current_rate);
                b[j] = (1 - 0.5 * pow(sigma_ * (j + 1), 2) * dt);
                c[j] = 0.25 * (j + 1) * dt * (pow(sigma_, 2) * (j + 1) + current_rate);
                d[j] = (1 + (current_rate + 0.5 * pow(sigma_ * (j + 1), 2)) * dt);
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
                V[i] = max(V[i], ds * (i + 1) - K_); // american call payoff 
                values_(N_ - n + 1, i + 1) = V[i];
            }

        }


        printMatrix();
    }

}