
/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Gabriel Espinoza Rojas
 * @Date  : 10.03.2018
 */




#ifndef ANPI_MULLER_HPP
#define ANPI_MULLER_HPP

#include "RootInterpolation.hpp"
#include "RootNewtonRaphson.hpp"

namespace anpi {

    using namespace std;

    const int MAX_ITERATIONS = 10000;
    template<typename T>
    T Muller(const boost::math::tools::polynomial<T>& poly,
                                              T a, const T eps) {

        T b = anpi::rootNewtonRaphson(poly, a, eps);
        T c = anpi::rootInterpolation(poly, a, b, eps);

        T res;

        for (int i = 0;;++i)     ///PREGUNTAR FOR
        {
            // Calculating various constants required
            // to calculate x3
            T f1 = poly.evaluate(a);
            T f2 = poly.evaluate(b);
            T f3 = poly.evaluate(c);
            T d1 = f1 - f3;
            T d2 = f2 - f3;
            T h1 = a - c;
            T h2 = b - c;
            T a0 = f3;
            T a1 = (((d2*h1*h1) - (d1*h2*h2))
                        / ((h1*h2) * (h1-h2)));
            T a2 = (((d1*h2) - (d2*h1))/((h1*h2) * (h1-h2)));
            T x = ((-2*a0) / (a1 + abs(std::sqrt(a1*a1-4*a0*a2))));
            T y = ((-2*a0) / (a1-abs(std::sqrt(a1*a1-4*a0*a2))));

            // Taking the root which is closer to x2
            if (x >= y)
                res = x + c;
            else
                res = y + c;

            // checking for resemblance of x3 with x2 till
            // two decimal places
            T m = res*T(100);
            float n = T*(100);
            m = floor(m);
            n = floor(n);
            if (m == n)
                break;
            a = b;
            b = c;
            c = res;
            if (i > MAX_ITERATIONS)
            {
                return std::numeric_limits<T>::quiet_NaN();
            }
        }
        if (i <= MAX_ITERATIONS){
            return res;
        }

        return std::numeric_limits<T>::quiet_NaN();
    }

}
#endif