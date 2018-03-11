/**
 * Copyright (C) 2017
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: Gabriel Espinoza
 * @Date  : 10.03.2018
 */

#include <boost/test/unit_test.hpp>

#include "Deflate.hpp"
#include "Deflate2.hpp"
#include "Muller.hpp"
#include "Laguerre.hpp"

#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>
#include <iostream>

#include <functional>

#include <cmath>

namespace anpi {
    namespace test {

        /// Square of a number
        template<typename T>
        inline T sqr(const T x) { return x*x; }

        /// First testing function for roots (x-1)(x²+ 9)(x+2)
        /// x^4 + x^3 + 7x^2 + 9x - 18
        template<typename T>
        T t1(const T x)  { return (x-1)*((x*x)+ 9)*(x+2); } /// raices reales y complejas x=1,\:x=3i,\:x=-3i,\:x=-2
        /// Second testing function for roots (x+1)(x²-9)(x+4)
        template<typename T>
        T t2(const T x)  { return (x+1)*((x*x)-9)*(x+2); }  /// 4 raices reales x=3,\:x=-1,\:x=-4,\:x=-3
        /// Third testing function for roots (x-2)(x²+16)(x+9)
        template<typename T>
        T t3(const T x)  { return (x-1)*((x*x)+ 16)*(x+9); } /// raices reales y complejas x=2,\:x=4i,\:x=-4i,\:x=-9



        enum TestIntervalMode {
            TestInterval,
            DoNotTestInterval
        };

        /// Test the given closed root finder
        template<typename T>
        void rootTest(const std::function<T(const std::function<T(T)>&,
                      T,
                      T,
                      const T)>& solver,
        const TestIntervalMode testInterval=TestInterval) {

        T eps=static_cast<T>(0.001);

        /* catch wrong interval */
        if (testInterval==TestInterval) {
        try {
        solver(t1<T>,T(2),T(0),eps);
        BOOST_CHECK(false && "solver should catch inverted interval");
    } catch(Exception exc ) {
    BOOST_CHECK(true && "successfully catched");
}

try {
solver(t3<T>,T(1),T(2),eps);
BOOST_CHECK(false && "solver should catch unenclosed root");
} catch(Exception exc ) {
BOOST_CHECK(true && "successfully catched");
}
}

for (T eps=T(1)/T(10); eps>static_cast<T>(1.0e-7); eps/=T(10)) {
T sol = solver(t1<T>, T(0), T(2), eps);
BOOST_CHECK(std::abs(t1<T>(sol)) < eps);
sol = solver(t2<T>, T(0), T(2), eps);
BOOST_CHECK(std::abs(t2<T>(sol))<eps);
sol = solver(t3<T>,T(0),T(0.5),eps);
BOOST_CHECK(std::abs(t3<T>(sol))<eps);
}
}

/// Test the given open root finder
template<typename T>
void rootTest(const std::function<T(const std::function<T(T)>&,
                                    T,
                                    const T)>& solver) {

    for (T eps=T(1)/T(10); eps>static_cast<T>(1.0e-7); eps/=T(10)) {
        T sol = solver(t1<T>,T(0),eps);
        BOOST_CHECK(std::abs(t1<T>(sol))<eps);
        sol = solver(t2<T>,T(2),eps);
        BOOST_CHECK(std::abs(t2<T>(sol))<eps);
        sol = solver(t3<T>,T(0),eps);
        BOOST_CHECK(std::abs(t3<T>(sol))<eps);
    }
}
} // test
}  // anpi

BOOST_AUTO_TEST_SUITE( RootFinder )

"Deflate.hpp"
"Deflate2.hpp"
 "Muller.hpp"
 "Laguerre.hpp"

BOOST_AUTO_TEST_CASE(Deflate)
        {

        }

BOOST_AUTO_TEST_CASE(Deflate2)
        {

        }

BOOST_AUTO_TEST_CASE(Muller)
        {

        }

BOOST_AUTO_TEST_CASE(Laguerre)
        {
            anpi::test::rootTest<float>(anpi::laguer);
            anpi::test::rootTest<double>(anpi::laguer);
        }


BOOST_AUTO_TEST_SUITE_END()