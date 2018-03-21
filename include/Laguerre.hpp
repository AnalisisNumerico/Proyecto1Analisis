/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Gabriel Espinoza
 * @Date  : 10.03.2018
 */

#ifndef ANPI_ROOT_SECANT_HPP
#define ANPI_ROOT_SECANT_HPP

///numerical recipes
#define MR 8
#define MT 10
#define MAXIT (MT*MR)



#include <boost/test/unit_test.hpp>
#include <boost/math/tools/polynomial.hpp>

namespace anpi {


    template<typename T>
    std::complex<T> laguer(const boost::math::tools::polynomial<T>& a, std::complex<T> x, T eps) {
        int m = a.degree();

        T abx, abp, abm, err;
        std::complex<T> dx, x1, b, d, f, g, h, sq, gp, gm, g2;
        static float frac[MR + 1] = {0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0};
        for (int iter = 1; iter<= MAXIT; iter++) {
            b = a[m];
            err = std::abs(b);
            d = f = std::complex<T>(0.0, 0.0);
            abx = std::abs(x);
            for (int j = m - 1; j>=0;j--) {
                f = ((x* f)+ d);
                d = ((x* d)+ b);
                b = ((x* b)+ a[j]);
                err = std::abs(b) + abx * err;
            }
            err *= eps;
            if (std::abs(b) <= err) {
                return x;           //We are on the root
            }
            g = (d/b);
            g2 = (g*g);
            h = (g2-(T(2.0)*(f/b)));
            sq = std::sqrt((float) (m - 1) * ((float) m* h)-g2);
            gp = g+sq;
            gm = g-sq;
            abp = std::abs(gp);
            abm = std::abs(gm);
            if (abp<abm)
                gp = gm;
            dx = ((std::fmax(abp, abm) > 0.0 ? (std::complex<T>((float) m, 0.0) / gp)
                                        : ((T(1) + abx) * std::complex<T>(cos((float) iter), sin((float) iter)))));
            x1 = x-dx;
            if (x.real() == x1.real() && x.imag() == x1.imag()) return x;    //converge
            if (iter % MT) x = x1;
            else x =(x-(frac[iter / MT]* dx));
        }
        return std::numeric_limits<T>::quiet_NaN();
    }
}

#endif
