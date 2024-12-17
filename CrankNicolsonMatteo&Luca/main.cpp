#include "Option.h"
#include "Diffmethod.h"
#include "Utils.h"  
#include "Matrix.h"
#include <iostream>
#include <iomanip> // For formatting

int main() {
    try {
        // Load option parameters
        m2::Option opt;
        opt.openInput("option_data.txt");

        opt.printParameters();

        if (opt.getAmerican())
        {
            m2::American ame(opt);

            if (opt.getCallPut())
            {
                ame.priceCall();
            }

            else { ame.pricePut();}

            ame.calculateDelta();

            std::cout << "Crank Nicolson price:" << std::endl;
            std::cout << ame.getPrice() << std::endl;

            m2::blackScholesPrice(opt.getCallPut(), opt.getS0(), opt.getK(), opt.getT(), m2::computeAverageRate(opt.getRates(), opt.getT()), opt.getSigma());

            ame.calculateDelta();
            ame.calculateGamma();
            ame.calculateTheta();
            ame.calculateVega();
            ame.calculateRho();

            std::cout << "Delta Crank Nicolson: " << ame.getDelta() << std::endl;
            std::cout << "Gamma Crank Nicolson: " << ame.getGamma() << std::endl;
            std::cout << "Theta Crank Nicolson: " << ame.getTheta() << std::endl;
            std::cout << "Vega Crank Nicolson: " << ame.getVega() << std::endl;
            std::cout << "Rho Crank Nicolson: " << ame.getRho() << std::endl;

            m2::writeOutputTxt(ame.getPrice(), ame.getDelta(), ame.getGamma(), ame.getTheta(), ame.getVega(), ame.getRho(), ame.getT0prices(), ame.getDeltaGraph(),ame.getBound());
        }

        else
        {
            m2::European eur(opt);

            if (opt.getCallPut())
            {
                eur.priceCall();
            }

            else { eur.pricePut(); }

            std::cout << "Crank Nicolson price:" << std::endl;
            std::cout << eur.getPrice() << std::endl;

            
            m2::blackScholesPrice(opt.getCallPut(), opt.getS0(), opt.getK(), opt.getT(), m2::computeAverageRate(opt.getRates(), opt.getT()), opt.getSigma());

            eur.calculateDelta();
            eur.calculateGamma();
            eur.calculateTheta();
            eur.calculateVega();
            eur.calculateRho();

            std::cout << "Delta Crank Nicolson: " << eur.getDelta() << std::endl;
            std::cout << "Gamma Crank Nicolson: " << eur.getGamma() << std::endl;
            std::cout << "Theta Crank Nicolson: " << eur.getTheta() << std::endl;
            std::cout << "Vega Crank Nicolson: " << eur.getVega() << std::endl;
            std::cout << "Rho Crank Nicolson: " << eur.getRho() << std::endl;

            m2::writeOutputTxt(eur.getPrice(), eur.getDelta(), eur.getGamma(), eur.getTheta(), eur.getVega(), eur.getRho(), eur.getT0prices(), eur.getDeltaGraph(), eur.getBound());

        }

     
        
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
