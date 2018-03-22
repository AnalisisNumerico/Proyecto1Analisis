/**
* Copyright (C) 2018
* Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
*
* This file is part of the CE3102 Numerical Analysis lecture at TEC
*
* @Author: JP
* @Date  : 10.03.2018
*/

#include <cstdlib>
#include <iostream>

#include "Deflate.hpp"
#include "Deflate2.hpp"

#include <boost/math/tools/polynomial.hpp>


void showErrorMessage() {
    std::cerr << "error"
              << std::endl;
}

/** main params defined as
 *  -m <l/m> -t <s/d> -e <epsilon> -r <root real part> -i <root imaginary part> -p <polynome>
 *
 *  -m selects the algortihm to be used to obtain the roots, default is muller
 *      m muller
 *      l laguerre
 *
 *  -t indicates the precision to be used, default is simple
 *      s simple
 *      d double
 *
 *  -e epsilon defines the precision threshold
 *
 *  -r the real part of the given root
 *
 *  -i the imaginary part of the given root
 *
 *  -p coeficients of the polynome, max size = 9
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char* argv[]) {

    if(argc > 13 && argc <= 23) {
        std::string metodo = std::string(argv[2]);
        std::string precision = std::string(argv[4]);
        if(!(precision.compare("d"))) {
            std::cout << "precision: double" << std::endl;
            double eps = atof(argv[6]);
            std::cout << "epsilon: " << eps << std::endl;
            double realRoot = atof(argv[8]);
            std::cout << "root real part: " << realRoot << std::endl;
            double imaginaryRoot = atof(argv[10]);
            std::cout << "root imaginary part: " << imaginaryRoot << std::endl;

            std::array<double,10> polyArray;
            for(int i = 12; i < argc; i++){
                polyArray[i-12] = atof(argv[i]);
            }
            boost::math::tools::polynomial<double>(polyArray.begin(),polyArray.size());
            std::cout << std::endl;

            for(int i = 0; i < argc-12; ++i)
                std::cout << polyArray[i] << ' ';

            if(!(metodo.compare("l"))){
                std::cout << "obtaining roots with laguerre algorithm..." << std::endl;
            }

            else {
                std::cout << "obtaining roots with muller algorithm..." << std::endl;
            }

        }
        else {
            std::cout << "used precision: simple" << std::endl;
            float eps = atof(argv[6]);
            std::cout << "epsilon: " << eps << std::endl;
            float realRoot = atof(argv[8]);
            std::cout << "root real part: " << realRoot << std::endl;
            float imaginaryRoot = atof(argv[10]);
            std::cout << "root imaginary part: " << imaginaryRoot << std::endl;

            std::array<float,10> polyArray;
            for(int i = 12; i < argc; i++){
                polyArray[i-12] = atof(argv[i]);
            }
            boost::math::tools::polynomial<float>(polyArray.begin(),polyArray.size());
            std::cout << std::endl;

            for(int i = 0; i < argc-12; ++i)
                std::cout << polyArray[i] << ' ';

            if(!(metodo.compare("l"))){
                std::cout << "obtaining roots with laguerre algorithm..." << std::endl;
            }

            else {
                std::cout << "obtaining roots with muller algorithm..." << std::endl;
            }
        }
    }
    else {
        showErrorMessage();
        return 1;
    }

    return 0;
}