#include <boost/math/tools/polynomial.hpp>

#ifndef ANPI_DEFLACION2
#define ANPI_DEFLACION2

namespace anpi {

    template<typename T>
    boost::math::tools::polynomial<T> deflate2(const boost::math::tools::polynomial<T>& poly,
                                               const std::complex<T>& root,
                                               boost::math::tools::polynomial<T>& residuo) { //lo estoy forzando a funcionar para complejos

        boost::math::tools::polynomial<double> q{{}};
        residuo = poly;
        q = poly;

        int n = poly.degree();

        for(int j = 0; j <= n; j++) {
            q[j] = T(0);
        }

        T raiz = abs(root);
        raiz *= raiz;

        for(int k = n - 2; k >= 0; k--) {
            q[k] = residuo[k + 2];
            residuo[k] -= q[k] * raiz; //talvez se quite
        }
        for(int i = 2; i <= n; i++) { //limpia el residuo
            residuo[i] = T(0);
        }

        return q;

    }
}

#endif