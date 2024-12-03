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

            std::cout << "B-S price:" << std::endl;
            std::cout << m2::blackScholesPrice(opt.getCallPut(), opt.getS0(), opt.getK(), opt.getT(), m2::computeAverageRate(opt.getRates(), opt.getT()), opt.getSigma()) << std::endl;

            ame.calculateDelta();

            std::cout << "Delta: " << ame.getDelta() << std::endl;
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

            std::cout << "B-S price:" << std::endl;
            std::cout << m2::blackScholesPrice(opt.getCallPut(), opt.getS0(), opt.getK(), opt.getT(), m2::computeAverageRate(opt.getRates(), opt.getT()), opt.getSigma()) << std::endl;

            eur.calculateDelta();

            std::cout << "Delta: " << eur.getDelta() << std::endl;

        }

     
        
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
