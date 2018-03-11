
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
#include "Laguerre.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/math/tools/polynomial.hpp>

int main() {

// Put your main code in here

    boost::array<double, 4> const lel = {{10, -6, -4, 3}};
    boost::math::tools::polynomial<double> const a(lel.begin(), (unsigned int) lel.end());

    std::cout << anpi::laguer(a, 0) << std::endl;


    return EXIT_FAILURE;
}