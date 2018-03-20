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
     *  using deflacion
     */
    BOOST_AUTO_TEST_CASE(Deflate_1){
        polinomial::polynomial<double> const poly{{-12,5,3}};
        polinomial::polynomial<double> r{{}};
        polinomial::polynomial<double> const q = {-4,3,0};
        std::array<double, 3> residueArray = {};
        residueArray.fill(0);
        polinomial::polynomial<double> const residue(residueArray.begin(), 2);
        double const root = -3;
        polinomial::polynomial<double> quotient = anpi::deflate<double>(poly,root,r);
        BOOST_CHECK(quotient == q);
        BOOST_CHECK(r == residue);
        std::cout << q << quotient << std::endl;
    }

    /** test for x³-2x²-4 divided by x-3
     *  using deflacion
     */
    BOOST_AUTO_TEST_CASE(Deflate_2){
        polinomial::polynomial<double> const poly{{-4,0,-2,1}};
        polinomial::polynomial<double> r{{}};
        polinomial::polynomial<double> const q{{3,1,1,0}};
        polinomial::polynomial<double> residue{{5,0,0,0}};
        double const root = 3;
        BOOST_CHECK((anpi::deflate<double>(poly,root,r)) == q);
        BOOST_CHECK(r == residue);
    }

    /** test for x⁵+5x⁴+9x³+11x²+12x+13 divided by x+2
     *  using deflacion
     */
    BOOST_AUTO_TEST_CASE(Deflate_3){
        polinomial::polynomial<double> const poly{{13,12,11,9,5,1}};
        polinomial::polynomial<double> r{{}};
        polinomial::polynomial<double> const q{{2,5,3,3,1,0}};
        polinomial::polynomial<double> residue{{9,0,0,0,0,0}};
        double const root = -2;
        BOOST_CHECK((anpi::deflate<double>(poly,root,r)) == q);
        BOOST_CHECK(r == residue);
    }

    /** test for x-1 divided by x+1
     *  using deflacion
     */
    BOOST_AUTO_TEST_CASE(Deflate_4){
        polinomial::polynomial<double> const poly{{1,1}};
        polinomial::polynomial<double> r{{}};
        polinomial::polynomial<double> const q{{1,0}};
        polinomial::polynomial<double> residue{{2,0}};
        double const root = 1;
        BOOST_CHECK((anpi::deflate<double>(poly,root,r)) == q);
        BOOST_CHECK(r == residue);
    }

BOOST_AUTO_TEST_SUITE_END()