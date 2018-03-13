/**
 * Copyright (C) 2017
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: JP
 * @Date  : 10.03.2018


#include <boost/test/unit_test.hpp>
#include <boost/math/tools/polynomial.hpp>

#include "Deflate.hpp"

#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>
#include <iostream>

#include <functional>

#include <cmath>

namespace anpi {
    namespace test {

        void deflacionTest1() {
            boost::math::tools::polynomial<double> const poly{{-12,5,3}};
            boost::math::tools::polynomial<double> r{{}};
            boost::math::tools::polynomial<double> const q{{-4,3,0}};
            double const root = -3;
            BOOST_CHECK((anpi::deflate<double>(poly,root,r)) == q);
        }

        void deflacionTest2() {
            boost::math::tools::polynomial<double> const poly{{double(-12), double(5), double(3)}};
            boost::math::tools::polynomial<double> r{{}};
            boost::math::tools::polynomial<double> const q{{double(-4), double(3) ,double(0)}};
            double const root = -3;
            BOOST_CHECK((anpi::deflate<double>(poly,root,r)) == q);
        }

        void deflacionTest3() {
            boost::math::tools::polynomial<double> const poly{{double(-12), double(5), double(3)}};
            boost::math::tools::polynomial<double> r{{}};
            boost::math::tools::polynomial<double> const q{{double(-4), double(3) ,double(0)}};
            double const root = -3;
            BOOST_CHECK((anpi::deflate<double>(poly,root,r)) == q);
        }

        void deflacionTest4() {
            boost::math::tools::polynomial<double> const poly{{double(-12), double(5), double(3)}};
            boost::math::tools::polynomial<double> r{{}};
            boost::math::tools::polynomial<double> const q{{double(-4), double(3) ,double(0)}};
            double const root = -3;
            BOOST_CHECK((anpi::deflate<double>(poly,root,r)) == q);
        }

    } // test
}  // anpi

BOOST_AUTO_TEST_SUITE( RootFinder )


    BOOST_AUTO_TEST_CASE(Deflate_1){
        anpi::test::deflacionTest1();
    }/*

    BOOST_AUTO_TEST_CASE(Deflate_2){
        anpi::test::deflacionTest2();
    }

    BOOST_AUTO_TEST_CASE(Deflate_3){
        anpi::test::deflacionTest3();
    }

    BOOST_AUTO_TEST_CASE(Deflate_4){
        anpi::test::deflacionTest4();
    }

BOOST_AUTO_TEST_SUITE_END()*/