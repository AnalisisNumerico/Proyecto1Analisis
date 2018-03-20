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

int main() {


    for(int j = 0; j < 1; j++) {
        boost::math::tools::polynomial<double> const poly{{-12, 5, 3}};
        boost::math::tools::polynomial<double> const q{{-4, 3, 0}};
        boost::math::tools::polynomial<double> r{{}};
        double const root = -3;

        std::cout << anpi::deflate<double>(poly,root,r);
        std::cout << r << std::endl;

        boost::math::tools::polynomial<double> const poly2{{1, 1, 1, 2}};
        std::complex<double> i;
        i = -1;
        std::complex<double> const root2 = sqrt(i);
        std::cout << anpi::deflate2(poly2,root2,r);
        std::cout << r;
        std::cout << anpi::deflate2(poly2,root2,r) << std::endl;

    }

    return 0;
}