#include "Option.h"
#include "CrankNicolson.h"
#include "utils.h"  
#include "Matrix.h"
#include <iostream>
#include <iomanip> // For formatting

int main() {
    try {
        // Load option parameters
        m2::Option opt;
        opt.openInput("C:/Users/matte/Documents/option_data.txt");

        opt.printParameters();

    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
