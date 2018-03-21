/**
 * Copyright (C) 2017
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: JP
 * @Date  : 10.03.2018*/


#include <boost/test/unit_test.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <iostream>
#include "Deflate.hpp"

namespace polinomial = boost::math::tools;

namespace anpi {
    namespace test {


    } // test
}  // anpi

BOOST_AUTO_TEST_SUITE( RootFinder )

    /** test for 3x²+5x-12 divided by x+3
     *  equals 3x-4 with a residue of 0
     *  using deflacion
     */
    BOOST_AUTO_TEST_CASE(Deflate_1){

        polinomial::polynomial<double> const poly{{-12,5,3}};
        double const root = -3;
        polinomial::polynomial<double> r{{}};
        polinomial::polynomial<double> const quotient = anpi::deflate<double>(poly,root,r);
        BOOST_CHECK(quotient[0] == -4);
        BOOST_CHECK(quotient[1] ==  3);
        BOOST_CHECK(quotient[2] ==  0);
        BOOST_CHECK(r[0] == 0);
    }

    /** test for x³-2x²-4 divided by x-3
     *  equals x²+x+3 with a residue of 5
     *  using deflacion
     */
       BOOST_AUTO_TEST_CASE(Deflate_2){
        polinomial::polynomial<double> const poly{{-4,0,-2,1}};
        double const root = 3;
        polinomial::polynomial<double> r{{}};
        polinomial::polynomial<double> quotient = anpi::deflate<double>(poly,root,r);
        BOOST_CHECK(quotient[0] == 3);
        BOOST_CHECK(quotient[1] == 1);
        BOOST_CHECK(quotient[2] == 1);
        BOOST_CHECK(quotient[3] == 0);
        BOOST_CHECK(r[0] == 5);
    }

    /** test for x⁵+5x⁴+9x³+11x²+12x+13 divided by x+2
     *  equals x⁴+3x³+3x²+5x+2 with a residue of 9
     *  using deflacion
     */
        BOOST_AUTO_TEST_CASE(Deflate_3){
        polinomial::polynomial<double> const poly{{13,12,11,9,5,1}};
        double const root = -2;
        polinomial::polynomial<double> r{{}};
        polinomial::polynomial<double> quotient = anpi::deflate<double>(poly,root,r);
        BOOST_CHECK(quotient[0] == 2);
        BOOST_CHECK(quotient[1] == 5);
        BOOST_CHECK(quotient[2] == 3);
        BOOST_CHECK(quotient[3] == 3);
        BOOST_CHECK(quotient[4] == 1);
        BOOST_CHECK(quotient[5] == 0);
        BOOST_CHECK(r[0] == 9);
    }

    /** test for x-1 divided by x+1
     *  equals 1 with a residue of 2
     *  using deflacion
     */
    BOOST_AUTO_TEST_CASE(Deflate_4){
        polinomial::polynomial<double> const poly{{1,1}};
        double const root = 1;
        polinomial::polynomial<double> r{{}};
        polinomial::polynomial<double> const q{{1,0}};
        polinomial::polynomial<double> residue{{2,0}};
        polinomial::polynomial<double> quotient = anpi::deflate<double>(poly,root,r);
        BOOST_CHECK(quotient[0] == 1);
        BOOST_CHECK(quotient[1] == 0);
        BOOST_CHECK(r[0] == 2);
        
    }

BOOST_AUTO_TEST_SUITE_END()