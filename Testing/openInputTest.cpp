#include <gtest/gtest.h>
#include "../CrankNicolsonMatteo&Luca/Option.h"
#include "../CrankNicolsonMatteo&Luca/Option.cpp"
#include <vector>
#include <utility>
#include <string>


// Test case for openInput
TEST(OptionTest, OpenInputValidFile) {
    // Arrange
    Option opt;

    // Act
    ASSERT_NO_THROW(opt.openInput("../test_parameters.txt"));

    // Assert
    // Check parameters after reading the file
    std::ostringstream output;
    opt.printParameters(); // Print the parameters for verification (can be removed)
    ASSERT_EQ(opt.getCallPut(), 0);  // Check correct call_put_
}

