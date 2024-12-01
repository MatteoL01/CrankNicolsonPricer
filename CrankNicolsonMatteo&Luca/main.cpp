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
        opt.openInput("C:/Users/matte/Documents/option_data.txt");

        opt.printParameters();

        m2::Matrix result;

        m2::American ame(opt);

        if (opt.getCallPut())
        {
            ame.priceCall();
        }

        else { ame.pricePut();}

        

    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
