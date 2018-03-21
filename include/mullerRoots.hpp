
/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Gabriel Espinoza Rojas
 * @Date  : 10.03.2018
 */

#ifndef ANPI_MULLER_ROOTS_HPP
#define ANPI_MULLER_ROOTS_HPP
#include <boost/test/unit_test.hpp>
#include <boost/math/tools/polynomial.hpp>


namespace anpi {

    // C++ Program to find root of a function, f(x)
#include<bits/stdc++.h>

    template<typename T>
    void mullerRoots(const boost::math::tools::polynomial<T>& polyOriginal, T x, std::vector<T>& roots, int polish, const T eps){
        int m = poly.degree();
        boost::math::tools::polynomial<T> poly=polyOriginal;
        boost::math::tools::polynomial<T> residuo;
        T newRoot;

        int j,i;
        for (j=0;j<=m;j++){
            newRoot = anpi::Muller(poly,x,eps);
            if (polish){
                    newRoot = anpi::Muller(polyOriginal,newRoot,eps);
            }
            roots.push_back(newRoot);
            std::cout << newRoot << std::endl;
            poly=anpi::deflate(poly, newRoot, residuo);
        }


    }
}
#endif