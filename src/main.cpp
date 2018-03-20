
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

using std::string;
using std::exception;
using std::cout;
using std::abs;
using std::pair;

using namespace boost::math::tools; // for polynomial



int main() {
/*
    boost::math::tools::polynomial<double> const poly{{double(-12), double(5), double(3)}};

    boost::math::tools::polynomial<double> const q{{double(-4), double(3) ,double(0)}};
    double const root = -3;
    std::cout << anpi::deflate<double>(poly,root,r);
    std::cout << r;
*/

    boost::array<double, 4> const d3a = {{10, -6, -4, 3}};
    polynomial<double> const a(d3a.begin(), 4);
    std::cout << "f1 " <<     a.evaluate(1) <<std::endl;


    return 0;
}