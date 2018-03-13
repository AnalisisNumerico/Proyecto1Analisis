#include <boost/math/tools/polynomial.hpp>

#ifndef ANPI_DEFLACION
#define ANPI_DEFLACION

namespace anpi {

    /**
     *
     * @param poly
     * @param root
     * @param residuo
     * @return
     */
    template<typename T>
    boost::math::tools::polynomial<T> deflate(const boost::math::tools::polynomial<T>& poly,
                                              const T& root,
                                              boost::math::tools::polynomial<T>& residuo) {

        boost::math::tools::polynomial<double> q{{}};
        residuo = poly;
        q = poly;

        int n = poly.degree();

        for(int j = 0; j <= n; j++) {
            q[j] = T(0);
        }

        for(int k = n - 1; k >= 0; k--) {
            q[k] = residuo[k + 1];
            residuo[k] -= q[k] * -root;
        }

        for(int i = 1; i <= n; i++) { //limpia el residuo
            residuo[i] = T(0);
        }

        return q;
    }
}

#endif