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
            std::cout << "Delta Crank Nicolson: " << eur.getDelta() << std::endl;
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
            std::cout << "Delta Crank Nicolson: " << eur.getDelta() << std::endl;
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
            std::cout << "Delta Crank Nicolson: " << ame.getDelta() << std::endl;
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
            std::cout << "Delta Crank Nicolson: " << ame.getDelta() << std::endl;
            std::cout << "----------------------------------------------------------------" << std::endl;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
