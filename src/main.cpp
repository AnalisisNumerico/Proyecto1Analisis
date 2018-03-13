
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

#include <boost/test/unit_test.hpp>
#include <boost/math/tools/polynomial.hpp>

using std::string;
using std::exception;
using std::cout;
using std::abs;
using std::pair;

using namespace boost::math;
using namespace boost::math::tools; // for polynomial
using boost::lexical_cast;



int main() {
/*
    boost::math::tools::polynomial<double> const poly{{double(-12), double(5), double(3)}};

    boost::math::tools::polynomial<double> const q{{double(-4), double(3) ,double(0)}};
    double const root = -3;
    std::cout << anpi::deflate<double>(poly,root,r);
    std::cout << r;
*/

    boost::math::tools::polynomial<double> const poly2{{1, 1, 1, 2}};
    boost::math::tools::polynomial<double> r{{}};
    std::complex<double> i;
    i = -1;
    std::complex<double> const root = sqrt(i);
    std::cout << anpi::deflate2(poly2,root,r);
    std::cout << r;

    return 0;
}