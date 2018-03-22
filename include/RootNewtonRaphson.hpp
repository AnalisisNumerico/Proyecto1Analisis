/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 */

#include <cmath>
#include <limits>
#include <functional>

#ifndef ANPI_NEWTON_RAPHSON_HPP
#define ANPI_NEWTON_RAPHSON_HPP

namespace anpi {

    /** Returns the evaluated derivate of a given function
     * by means of forward aproximation
     *
     *
     * @param funct a functor of the form "T funct(T x)"
     * @param value
     * @return evalueted derivate aproximation
     */
    template<typename T>
    T primeraDerivada(const boost::math::tools::polynomial<float>& funct, T x,T eps) {
        const T h = std::abs(eps) / T(2);
        return ((funct.evaluate(x+h) - funct.evaluate(x))/h);
    }

    /*template<typename T>
    std::complex<T> primeraDerivada(const boost::math::tools::polynomial<std::complex<T>>& funct, std::complex<T> x,T eps) {
        const std::complex<T> h = std::abs(eps) / std::complex<T>(2,0);
        return ((funct.evaluate(x+h) - funct.evaluate(x))/h);
    }*/

    /**
     * Find the roots of the function funct looking by means of the
     * Newton-Raphson method
     *
     * @param funct a functor of the form "T funct(T x)"
     * @param xi initial root guess
     *
     * @return root found, or NaN if none could be found.
     *
     * @throws anpi::Exception if inteval is reversed or both extremes
     *         have same sign.
     */
    template<typename T>
    T rootNewtonRaphson(const boost::math::tools::polynomial<T>& funct,T xi,const T eps) {

        int const MAX_ITERATIONS = 20;

        T x = xi;
        T dx;

        for(int i = 0; i < MAX_ITERATIONS; i++) {
            dx = ((funct.evaluate(x))/(primeraDerivada(funct,x,eps)));
            x = x - dx;
            if(std::abs(dx) < eps) {
                return x;
            }
        }

        // Return NaN if no root was found
        return std::numeric_limits<T>::quiet_NaN();
    }

    /*template<typename T>
    std::complex<T> rootNewtonRaphson(const boost::math::tools::polynomial<std::complex<T>>& funct,T xi,const T eps) {

        int const MAX_ITERATIONS = 20;

        std::complex<T> x = xi;
        std::complex<T> dx;

        for(int i = 0; i < MAX_ITERATIONS; i++) {
            dx = ((funct.evaluate(x))/(primeraDerivada(funct,x,eps)));
            x = x - dx;
            if(std::abs(dx) < eps) {
                return x;
            }
        }

        // Return NaN if no root was found
        return std::numeric_limits<T>::quiet_NaN();
    }*/

}

#endif