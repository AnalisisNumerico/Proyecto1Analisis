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

#include "Deflate.hpp"

#include <iostream>

namespace anpi {
    namespace test {

        /** test for 3x²+5x-12 divided by x+3
         *  using deflacion
         */
        void deflacionTest1() {
            boost::math::tools::polynomial<double> const poly{{-12,5,3}};
            boost::math::tools::polynomial<double> r{{}};
            boost::math::tools::polynomial<double> const q{{-4,3,0}};
            boost::math::tools::polynomial<double> residue{{0,0,0}};
            double const root = -3;
            BOOST_CHECK((anpi::deflate<double>(poly,root,r)) == q);
            BOOST_CHECK(r == residue);

            std::cout << "test1 passed" << std::endl;
        }

        /** test for x³-2x²-4 divided by x-3
         *  using deflacion
         */
        void deflacionTest2() {
            boost::math::tools::polynomial<double> const poly{{-4,0,-2,1}};
            boost::math::tools::polynomial<double> r{{}};
            boost::math::tools::polynomial<double> const q{{3,1,1,0}};
            boost::math::tools::polynomial<double> residue{{5,0,0,0}};
            double const root = 3;
            BOOST_CHECK((anpi::deflate<double>(poly,root,r)) == q);
            BOOST_CHECK(r == residue);
        }

        /** test for x⁵+5x⁴+9x³+11x²+12x+13 divided by x+2
         *  using deflacion
         */
        void deflacionTest3() {
            boost::math::tools::polynomial<double> const poly{{13,12,11,9,5,1}};
            boost::math::tools::polynomial<double> r{{}};
            boost::math::tools::polynomial<double> const q{{2,5,3,3,1,0}};
            boost::math::tools::polynomial<double> residue{{9,0,0,0,0,0}};
            double const root = -2;
            BOOST_CHECK((anpi::deflate<double>(poly,root,r)) == q);
            BOOST_CHECK(r == residue);
        }

        /** test for x-1 divided by x+1
         *  using deflacion
         */
        void deflacionTest4() {
            boost::math::tools::polynomial<double> const poly{{1,1}};
            boost::math::tools::polynomial<double> r{{}};
            boost::math::tools::polynomial<double> const q{{1,0}};
            boost::math::tools::polynomial<double> residue{{2,0}};
            double const root = 1;
            BOOST_CHECK((anpi::deflate<double>(poly,root,r)) == q);
            BOOST_CHECK(r == residue);
        }

    } // test
}  // anpi

BOOST_AUTO_TEST_SUITE( RootFinder )

    BOOST_AUTO_TEST_CASE(Deflate_1){
        std::cout << "test1 passed" << std::endl;
        anpi::test::deflacionTest1();
    }

    BOOST_AUTO_TEST_CASE(Deflate_2){
        anpi::test::deflacionTest2();
    }

    BOOST_AUTO_TEST_CASE(Deflate_3){
        anpi::test::deflacionTest3();
    }

    BOOST_AUTO_TEST_CASE(Deflate_4){
        anpi::test::deflacionTest4();
    }

BOOST_AUTO_TEST_SUITE_END()