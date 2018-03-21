/**
* Copyright (C) 2018
* Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
*
* This file is part of the numerical analysis lecture CE3102 at TEC
*
* @Author: Gabriel Espinoza Rojas
* @Date : 10.03.2018
*/




#ifndef ANPI_MULLER_HPP
#define ANPI_MULLER_HPP

#include "RootInterpolation.hpp"
#include "RootNewtonRaphson.hpp"
#include <iostream>

namespace anpi {

    const int MAX_ITERATIONS = 10000;
    template<typename T>
    T Muller(const boost::math::tools::polynomial<T>& poly,
             T a, const T eps) {
        T b;
        T c;
        b = anpi::rootNewtonRaphson(poly, a, eps);
        if(a<b){
            c = anpi::rootInterpolation(poly, a, b, eps);
        }else{
            c = anpi::rootInterpolation(poly, b, a, eps);
        }

        //std::cout << "a " <<a<< std::endl;
        //std::cout << "b " << b<<std::endl;
        //std::cout << "c " << c<<std::endl;

        //a = T(0);
        //b = T(-0.505356);
        //c = T(-0.5);




        T res;
        int i;

        for (i = 0;;++i){
            T f1 = poly.evaluate(a); // 1
            T f2 = poly.evaluate(b); //1
            T f3 = poly.evaluate(c); //1
            T d1 = f1 - f3; //0
            T d2 = f2 - f3; //0
            T h1 = a - c; //0.50006  // xi - xi-1
            T h2 = b - c; //0,000075 // xi-1 - xi-2
            T a0 = f3; //1
            T a1 = (((d2*h1*h1) - (d1*h2*h2))
                    / ((h1*h2) * (h1-h2))); //0
            T a2 = (((d1*h2) - (d2*h1))/((h1*h2) * (h1-h2))); //0
            T x = ((-T(2)*a0) / (a1 + abs(std::sqrt(a1*a1-T(4)*a0*a2)))); //
            //std::cout << "x " << x<<std::endl;
            T y = ((-T(2)*a0) / (a1-abs(std::sqrt(a1*a1-T(4)*a0*a2)))); //
            //std::cout << "y " << y<<std::endl;

            if (x >= y)
                res = x + c;
            else
                res = y + c;

            T m = res*T(100);
            float n = c*T(100);
            m = floor(m);
            n = floor(n);
            if (m == n)
                break;
            a = b;
            b = c;
            c = res;
            if (i > MAX_ITERATIONS){
                std::cout << "res " << res<<std::endl;
                return std::numeric_limits<T>::quiet_NaN();
            }
        }
        if (i <= MAX_ITERATIONS){
            return res;
        }
        std::cout << "res " << res<<std::endl;
        return std::numeric_limits<T>::quiet_NaN();
    }

    template<typename T>
    void mullerRoots(const boost::math::tools::polynomial<T>& polyOriginal, T x, std::vector<T>& roots, int polish, const T eps){
        int m = polyOriginal.degree();
        boost::math::tools::polynomial<T> poly=polyOriginal;
        boost::math::tools::polynomial<T> residuo;
        T newRoot;

        int j;
        for (j=0;j<=m;j++){
            newRoot = anpi::Muller(poly,x,eps);
            if (polish){
                newRoot = anpi::Muller(polyOriginal,newRoot,eps);
            }
            roots.push_back(newRoot);
            std::cout << newRoot << std::endl;

            poly=anpi::deflate(poly, newRoot, residuo);
        }


    }

}
#endif