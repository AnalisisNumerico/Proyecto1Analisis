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
#include <boost/math/tools/polynomial.hpp>

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


        ///x^4-5x^2+5x^3-45x+4x^2-36
        /// First testing function for roots (x-1)(x²+ 9)(x+2)
        /// x^4 + x^3 + 7x^2 + 9x - 18

        template<typename T>
        boost::math::tools::polynomial<T> t1(){
            boost::array<T, 5> const a = {{-18.0, 9.0,7.0, 1.0,1.0}};

            boost::math::tools::polynomial<T> const t(a.begin(), 5); /// raices reales y complejas x=1,\:x=3i,\:x=-3i,\:x=-2
            return t;
        }


        template<typename T>
        /// Second testing function for roots (x+1)(x²-9)(x+4)
        ///x^4+3x^3-7x^2-27x-18
        boost::math::tools::polynomial<T> t2(){
            boost::array<T, 5> const a = {{-18.0,-27.0,-7.0, 3.0,1.0}};
            boost::math::tools::polynomial<T> const t(a.begin(), 5); /// 4 raices reales x=-1,\:x=-2,\:x=-3,\:x=3
            return t;
        }

        template<typename T>
        /// Third testing function for roots (x-2)(x²+16)(x+9)
        ///x^4+7x^3-2x^2+112x-288
        boost::math::tools::polynomial<T> t3(){
            boost::array<T, 5> const a = {{-288.0,112.0,-2.0, 7.0,1.0}};

            boost::math::tools::polynomial<T> const t(a.begin(), 5); /// raices reales y complejas x=2,\:x=4i,\:x=-4i,\:x=-9
            return t;
        }



        enum TestIntervalMode {         ///COMO USAR ESTO PARA EL PULIDO
            TestInterval,
            DoNotTestInterval
        };

        /// Test the given a root finder
        template<typename T>
        void rootTest() {
            //T eps=T(1)/T(10);
            //std::vector<T> rootst1;
            //anpi::mullerRoots(t1<T>(), T(0), rootst1,0,eps);

            T eps=T(1)/T(10);
            std::complex<T> sol = anpi::laguer(t1<T>(), std::complex<T>(0.0,0.0), eps);
            std::cout << "Sol t1 "<< sol << std::endl;
            sol = anpi::laguer(t2<T>(), std::complex<T>(0.0,0.0), eps);
            std::cout << "Sol t2 "<< sol << std::endl;
            sol = anpi::laguer(t3<T>(), std::complex<T>(0.0,0.0), eps);
            std::cout << "Sol t3 "<< sol << std::endl;


            /*
            T eps=T(1)/T(10);
            T sol = anpi::Muller(t1<T>(), T(0), eps);
            std::cout << "Sol t1 "<< sol << std::endl;
            sol = anpi::Muller(t2<T>(), T(0), eps);
            std::cout << "Sol t2 "<< sol << std::endl;
            sol = anpi::Muller(t3<T>(), T(0), eps);
            std::cout << "Sol t3 "<< sol << std::endl;

            T valorTeoricot1 = T(1);
            T valorTeoricot2 = T(-1);
            T valorTeoricot3 = T(2);

            for (T eps=T(1)/T(10); eps>static_cast<T>(1.0e-7); eps/=T(10)) {

                T sol = anpi::Muller(t1<T>(), T(0), eps);
                //T sol = anpi::Muller(t1<T>(), std::complex<T>(0.0,0.0), eps);
                T error = valorTeoricot1 - sol;
                std::cout << "Sol t1 "<< sol << " eps "<< eps << std::endl;
                std::cout << "err t1 "<< error << " eps "<< eps << std::endl;
                BOOST_CHECK(error < eps);
                sol = anpi::Muller(t2<T>(), T(0), eps);
                error = valorTeoricot2 - sol;
                //BOOST_CHECK(error<eps);
                sol = anpi::Muller(t3<T>(), T(0), eps);
                //error = valorTeoricot3 - sol;
                BOOST_CHECK(error<eps);
            }*/

        }


        /*
       /// Test the given all the roots
       template<typename T>                                                            ///ES UNA STD::FUNCTION EL SOLVER?
       void rootsTest(const std::function<T(const boost::math::tools::polynomial<T>&,
                                            T, std::vector<T>&, int, const T)>& solver,
                      const TestIntervalMode testInterval=TestInterval) {
           std::vector<T> rootst1;
           std::vector<T> rootst2;
           std::vector<T> rootst3;
           std::vector<T> rootsSolt1 = solver(t1<T>, T(0),rootst1, 0,eps);
           std::vector<T> rootsSolt2 = solver(t2<T>, T(0),rootst2, 0,eps);
           std::vector<T> rootsSolt3 = solver(t3<T>, T(0),rootst3, 0,eps);
           for(int i=0; i<rootsSolt1.size(); i++) {
               for (T eps=T(1)/T(10); eps>static_cast<T>(1.0e-7); eps/=T(10)) {
                   BOOST_CHECK(std::abs(t1<T>(rootsSolt1[i])) < eps);
                   BOOST_CHECK(std::abs(t2<T>(rootsSolt1[i]))<eps);
                   BOOST_CHECK(std::abs(t3<T>(rootsSolt1[i]))<eps);
               }
           }
       }*/

    } // test
}  // anpi


BOOST_AUTO_TEST_SUITE( RootFinder )


    BOOST_AUTO_TEST_CASE(Deflate)
    {

    }

    BOOST_AUTO_TEST_CASE(Deflate2)
    {

    }

    BOOST_AUTO_TEST_CASE(Muller)
    {
        anpi::test::rootTest<float>();
        anpi::test::rootTest<double>();
    }

    BOOST_AUTO_TEST_CASE(Laguerre)
    {
//anpi::test::rootTest<float>(anpi::laguerre);
        //anpi::test::rootTest<double>(anpi::laguerre);
    }
    BOOST_AUTO_TEST_CASE(mullerRoots){
        //anpi::test::rootTest<float>(anpi::mullerRoots);
        //anpi::test::rootTest<double>(anpi::mullerRoots);
    }
    BOOST_AUTO_TEST_CASE(laguerreRoots){

        //anpi::test::rootTest<float>(anpi::laguerreRoots);
        //anpi::test::rootTest<double>(anpi::laguerreRoots);
    }

BOOST_AUTO_TEST_SUITE_END()