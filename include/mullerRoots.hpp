
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
#include "Muller.hpp"
#include "Deflate.hpp"



namespace anpi {

    // C++ Program to find root of a function, f(x)
#include<bits/stdc++.h>
    using namespace std;

    const int MAX_ITERATIONS = 10000;

    template<typename T>
    void mullerRoots(const boost::math::tools::polynomial<T>& polyOriginal, int x, const std::vector<T>& roots, int polish, const T eps) {
        m = poly.degree();
        std::complex x,b,c,ad[MAXM];
        boost::math::tools::polynomial<T> poly=polyOriginal;
        boost::math::tools::polynomial<T> residuo;
        for (int j=0;j<=m;j++){
            newRoot = anpi::Muller(poly,x,eps);
            roots.push_back(newRoot);
            poly=anpi::deflate(poly, x, residuo);
            //x=newRoot;
        }

        if (polish)
            for (j=1;j<=m;j++)
                anpi::Muller(poly,x,eps);
            for (j=2;j<=m;j++){
                for (i=j-1;i>=1;i--) {
                    if (roots[i].r <= x.r) break;
                    roots[i+1]=roots[i];
                }
                roots[i+1]=x;
            }
    }
}
#endif