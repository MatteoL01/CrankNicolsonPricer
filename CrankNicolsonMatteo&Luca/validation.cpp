#include "Option.h"
#include "Diffmethod.h"
#include "Utils.h"  
#include "Matrix.h"
#include <iostream>
#include <iomanip> // For formatting

int main() {
    try {
        // Test for the European Call
        {
            m2::Option opt;
            opt.openInput("test_Ecall.txt");

            opt.printParameters();

            m2::European eur(opt);
            eur.priceCall();
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

            std::cout << "----------------------------------------------------------------" << std::endl;
        }

        // Test for the European Put
        {
            m2::Option opt;
            opt.openInput("test_Eput.txt");

            opt.printParameters();

            m2::European eur(opt);
            eur.pricePut();
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

            std::cout << "----------------------------------------------------------------" << std::endl;
        }

        // Test for the American Call
        {
            m2::Option opt;
            opt.openInput("test_Acall.txt");

            opt.printParameters();

            m2::American ame(opt);
            ame.priceCall();
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

            std::cout << "----------------------------------------------------------------" << std::endl;
        }

        // Test for the American Put
        {
            m2::Option opt;
            opt.openInput("test_Aput.txt");

            opt.printParameters();

            m2::American ame(opt);
            ame.pricePut();
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

            std::cout << "----------------------------------------------------------------" << std::endl;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
