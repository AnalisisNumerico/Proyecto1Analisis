#include <boost/math/tools/polynomial.hpp>

#ifndef ANPI_DEFLACION2
#define ANPI_DEFLACION2

namespace polinomial = boost::math::tools;

namespace anpi {

    template<typename T>
    polinomial::polynomial<T> deflate2(const polinomial::polynomial<T>& poly,
                                       const std::complex<T>& root,
                                       polinomial::polynomial<T>& residuo) {

        polinomial::polynomial<double> q{{}};
        residuo = poly;
        q = poly;

        int n = poly.degree();

        for(int j = 0; j <= n; j++) {
            q[j] = T(0);
        }

        T raiz2 = 2 * root.real();

        T raiz1 = norm(root);

        for(int k = n - 2; k >= 0; k--) {
            q[k] = residuo[k + 2];
            residuo[k+1] -= q[k] * raiz2;
            residuo[k] -= q[k] * raiz1;
        }
        for(int i = 2; i <= n; i++) { //limpia el residuo
            residuo[i] = T(0);
        }

        return q;

    }
}

#endif