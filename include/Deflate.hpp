#include <boost/math/tools/polynomial.hpp>

#ifndef ANPI_DEFLACION
#define ANPI_DEFLACION

namespace polinomial = boost::math::tools;

namespace anpi {

    /**
     *
     * @param poly
     * @param root
     * @param residuo
     * @return
     */
    template<typename T>
    polinomial::polynomial<T> deflate(const polinomial::polynomial<T>& poly,
                                      const T& root,
                                      polinomial::polynomial<T>& residuo) {

/*
        residuo = poly;
        boost::math::tools::polynomial<T> q = {{}};
        q = poly;
        int n = poly.degree();
        for(int j = 0; j <= n; j++) {
            q[j] = T(0);
        }
        for(int k = n - 1; k >= 0; k--) {
            q[k] = residuo[k + 1];
            residuo[k] -= q[k] * -root;
        }
        for(int i = 1; i <= n; i++) {
            residuo[i] = T(0);
        }
        return q;
*/
        boost::array<T, 0> const a = {{}};
        boost::math::tools::polynomial<T> const q(a.begin(), 0);
        //polinomial::polynomial<T> q = {{}};

        q = poly;
        int n = poly.degree();
        T remanente = poly[n];
        q[n] = T(0);

        for(int i = n - 1; i >= 0; i--){
            T swap = q[i];
            q[i] = remanente;
            remanente = swap + remanente * root;
        }

        residuo = poly;
        residuo[0] = remanente;

        for(int i = 1; i <= n; i++) {
            residuo[i] = T(0);
        }

        return q;
    }
}

#endif