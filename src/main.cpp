/**
* Copyright (C) 2018
* Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
*
* This file is part of the CE3102 Numerical Analysis lecture at TEC
*
* @Author: Gabriel Espinoza Rojas
* @Date  : 10.03.2018
*/

#include <cstdlib>
#include <iostream>
#include <vector>
#include <complex>
#include <iomanip>

#include "Deflate.hpp"
#include "Deflate2.hpp"

#include <boost/math/tools/polynomial.hpp>


int main(int argc, char* argv[]) {
    // Check the number of parameters
    if (argc < 2) {
        // Tell the user how to run the program
        std::cerr << "Usage: " << argv[0] << " NAME" << std::endl;
        /* "Usage messages" are a conventional way of telling the user
         * how to run a program if they enter the command incorrectly.
         */
        return 1;
    }
    // Print the user's name:
    std::cout << argv[0] << "says hello, " << argv[1] << "!" << std::endl;
    return 0;
}