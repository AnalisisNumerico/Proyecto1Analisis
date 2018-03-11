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
#define EPSS 1.0e-7
#define MR 8
#define MT 10
#define MAXIT (MT*MR)

//Here EPSS is the estimated fractional roundoff error. We try to break (rare) limit cycles with MR
//different fractional values, once every MT steps, for MAXIT total allowed iterations.

#include <boost/test/unit_test.hpp>
#include <boost/math/tools/polynomial.hpp>

namespace anpi {

    /**
     * The following routine implements the Laguerre method to find one root of a given polynomial of degree m,
     * whose coefficients can be complex. As usual, the first coefficient a[0] is the constant term, while a[m]
     * is the coefficient of the highest power of x. The routine implements a simplified version of an elegant
     * stopping criterion due to Adams[5], which neatly balances the desire to achieve full machine
     * accuracy, on the one hand, with the danger of iterating forever in the presence of
     * roundoff  error,  on  the  other.
     *
     * @Author: Gabriel Espinoza
     * @Date  : 10.03.2018
     * @param m es el grado del polinomio
     */
    template<typename T>

    T laguer(const boost::math::tools::polynomial<T>& a, std::complex<T> *x) {
        int m = a.degree();
        float abx, abp, abm, err;
        std::complex<T> dx, x1, b, d, f, g, h, sq, gp, gm, g2;
        static float frac[MR + 1] = {0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0};
        for (int iter = 1; iter<= MAXIT; iter++) {
            b = a[m];
            err = Cabs(b);
            d = f = Complex(0.0, 0.0);
            abx = Cabs(*x);
            for (int j = m - 1; j>=0;j--) {
                f = Cadd(Cmul(*x, f), d);
                d = Cadd(Cmul(*x, d), b);
                b = Cadd(Cmul(*x, b), a[j]);
                err = Cabs(b) + abx * err;
            }
            err *= EPSS;
            if (Cabs(b) <= err) return x;       //We are on the root
            g = Cdiv(d, b);
            g2 = Cmul(g, g);
            h = Csub(g2, RCmul(2.0, Cdiv(f, b)));
            sq = Csqrt(RCmul((float) (m - 1), Csub(RCmul((float) m, h), g2)));
            gp = Cadd(g, sq);
            gm = Csub(g, sq);
            abp = Cabs(gp);
            abm = Cabs(gm);
            if (abp<abm)
                gp = gm;
            dx = ((FMAX(abp, abm) > 0.0 ? Cdiv(Complex((float) m, 0.0), gp)
                                        : RCmul(1 + abx, Complex(cos((float) iter), sin((float) iter)))));
            x1 = Csub(*x, dx);
            if (x->r == x1.r &&x->i == x1.i) return x;    //converge
            if (iter % MT) *x = x1;
            else *x = Csub(*x, RCmul(frac[iter / MT], dx));
        }
        nrerror("too many iterations in laguer");

        return std::numeric_limits<T>::quiet_NaN();
    }
}

#endif
