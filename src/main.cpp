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

int main(int argc, char* argv[]) {

    if(argc > 2 && strcmp(argv[1],"-p")) {

        int polyBegin = 2;
        int polyEnd = -1;
        std::string metodo = "x";
        std::string precision = "x";
        double eps = -1;

        for(int i = 2; i < argc; i++) {
            if(strcmp(argv[i],"-m")) {
                polyEnd = i;
                if (i + 5 < argc) {
                    metodo = std::string(argv[i + 1]);
                    precision = std::string(argv[i + 3]);
                    eps = atof(argv[i + 5]);
                }
            }
        }

        if(polyEnd == -1 || metodo.compare("x") ||
            precision.compare("x") || eps == -1) {
            showErrorMessage();
            return 1;
        }

        for(int i = polyBegin; i < polyEnd; i++){
        }

    }
    else {
        showErrorMessage();
        return 1;
    }

    return 0;
}