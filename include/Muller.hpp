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
    std::complex<T> Muller(const boost::math::tools::polynomial<T>& poly,
             std::complex<T>& xi, const T eps) {

        T a;
        T b;
        T c;

        if(xi.imag()==0){
            a=xi.real();
            b = anpi::rootNewtonRaphson(poly, a, eps);
            if(a<b){
                std::cout << "1a " <<a<< std::endl;
                std::cout << "1b " << b<<std::endl;
                c = anpi::rootInterpolation(poly, a, b, eps);
            }else{
                std::cout << "2a " <<a<< std::endl;
                std::cout << "2b " << b<<std::endl;
                c = anpi::rootInterpolation(poly, b, a, eps);
            }
        }else{


            std::cout << "a " << std::endl;
            //std::complex<T> b = (br,0);
            //std::complex<T> c = (cr,0);
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
            T f1, f2, f3, d1, d2, h1, h2, a0, a1, a2, x, y;
            f1 = poly.evaluate(a); // 1
            f2 = poly.evaluate(b); //1
            f3 = poly.evaluate(c); //1
            d1 = f1 - f3; //0
            d2 = f2 - f3; //0
            h1 = a - c; //0.50006  // xi - xi-1
            h2 = b - c; //0,000075 // xi-1 - xi-2
            a0 = f3; //1
            a1 = (((d2*h1*h1) - (d1*h2*h2))
                    / ((h1*h2) * (h1-h2))); //0
            a2 = (((d1*h2) - (d2*h1))/((h1*h2) * (h1-h2))); //0
            x = ((-T(2)*a0) / (a1 + abs(std::sqrt(a1*a1-T(4)*a0*a2)))); //
            //std::cout << "x " << x<<std::endl;
            y = ((-T(2)*a0) / (a1-abs(std::sqrt(a1*a1-T(4)*a0*a2)))); //
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
    std::complex<T> Muller(const boost::math::tools::polynomial<std::complex<T>>& poly,
                           std::complex<T>& a, const T eps) {


//        std::complex<T> b (a.real()+1, a.imag()+1);
//        std::complex<T> c (a.real()+2, a.imag()+2);

        std::complex<T> b (1, 1);
        std::complex<T> c (2, 2);


        //std::cout << "a " << a << std::endl;
        //std::cout << "b " << b <<std::endl;
        //std::cout << "c " << c <<std::endl;

        //a = T(0);
        //b = T(-0.505356);
        //c = T(-0.5);




        std::complex<T> res;
        int i;



        for (i = 0;;++i){
            std::complex<T> f1, f2, f3, d1, d2, h1, h2, a0, a1, a2;
            T x, y;

            f1 = poly.evaluate(a); // 1
            f2 = poly.evaluate(b); //1
            f3 = poly.evaluate(c); //1
            d1 = f1 - f3; //0
            d2 = f2 - f3; //0
            h1 = a - c; //0.50006  // xi - xi-1
            h2 = b - c; //0,000075 // xi-1 - xi-2
            a0 = f3; //1
            a1 = (((d2*h1*h1) - (d1*h2*h2))
                  / ((h1*h2) * (h1-h2))); //0
            a2 = (((d1*h2) - (d2*h1))/((h1*h2) * (h1-h2))); //0

            x = ((-T(2)*a0) / (a1 + abs(std::sqrt(a1*a1-T(4)*a0*a2)))); //
            y = ((-T(2)*a0) / (a1-abs(std::sqrt(a1*a1-T(4)*a0*a2)))); //
            std::cout << "x "<<std::endl;

            x = ((std::complex<T>(-2,0)*a0) / (a1 + abs(std::sqrt(a1*a1-std::complex<T>(4,0)*a0*a2)))); //
            std::cout << "x "<<std::endl;
            y = ((std::complex<T>(-2,0)*a0) / (a1-abs(std::sqrt(a1*a1-std::complex<T>(4,0)*a0*a2)))); //

            //std::cout << "y " << y<<std::endl;

            if (std::abs(x) >= std::abs(y))
                res = x + c;
            else
                res = y + c;

            T m = res.real()*T(100);
            float n = c.real()*T(100);
            m = floor(m);
            n = floor(n);
            if (m == n)
                break;
            a = b;
            b = c;
            c = res;
            if (i > MAX_ITERATIONS){
                std::cout << "res " << res<<std::endl;
                return res;
            }
        }
        if (i <= MAX_ITERATIONS){
            return res;
        }
        std::cout << "res " << res<<std::endl;
        return res;
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