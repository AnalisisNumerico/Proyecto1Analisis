/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Gabriel Espinoza
 * @Date  : 10.03.2018
 */

#ifndef ANPI_LAGUERRE_HPP
#define ANPI_LAGUERRE_HPP





#include <boost/test/unit_test.hpp>
#include <boost/math/tools/polynomial.hpp>

namespace anpi {


    template<typename T>
    std::complex<T> laguer(const boost::math::tools::polynomial<T>& a, std::complex<T>& x, T eps) {
        T const MR = 8;
        std::complex<T> const MT = 10;
        int const MAXIT = 80;
        int m = a.degree();


        T abx, abp, abm, err;
        std::complex<T> dx, x1, b, d, f, g, h, sq, gp, gm, g2;
        static float frac[10 + 1] = {0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0};
        for (int iter = 1; iter<= MAXIT; iter++) { //Loop over iterations up to allowed maximum.
            b = a[m];
            err = std::abs(b);
            d = f = std::complex<T>(0.0, 0.0);
            abx = std::abs(x);
            for (int j = m - 1; j >= 0; j--) { //Efficient computation of the polynomial and its first two derivatives. f stores P'' /2.
                f = ((x * f) + d);
                d = ((x * d) + b);
                b = ((x * b) + a[j]);
                err = std::abs(b) + abx * err;
            }
            err *= eps; //Estimate of roundoff error in evaluating polynomial.
            if (std::abs(b) <= err) {
                return x;           //We are on the root
            }
            g = (d / b);            //The generic case: use Laguerre’s formula.
            g2 = (g * g);
            h = (g2 - (T(2.0) * (f / b)));
            sq = sqrt(((T)(m-1))*(((T)(m)*h)-g2));
            gp = g + sq;
            gm = g - sq;
            abp = std::abs(gp);
            abm = std::abs(gm);

            if (abp < abm) {
            gp = gm;
            }

            dx=((fmax(abp,abm) > 0.0 ? ((std::complex<T>(m,0.0))/gp) : (1+abx)*(std::complex<T>(cos(iter),sin(iter)))));

            x1=x-dx;

            if(x == x1) return x;


            if (iter % (int) MT.real()){
                x = x1;
            }
            else {
                x = (x - (frac[iter / 10] * dx.real()));
            }
        }

        //Every so often we take a fractional step, to break any limit cycle (itself a rare occur-rence).
        return std::numeric_limits<T>::quiet_NaN();
    }
}

#endif
