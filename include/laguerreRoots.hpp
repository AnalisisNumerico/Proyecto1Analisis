
/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Gabriel Espinoza Rojas
 * @Date  : 10.03.2018
 */




#ifndef ANPI_LAGUERRE_ROOTS_HPP
#define ANPI_LAGUERRE_ROOTS_HPP
#include "Laguerre.hpp"


namespace anpi {

    // C++ Program to find root of a function, f(x)
#include<bits/stdc++.h>
    using namespace std;


    template<typename T>
    void laguerreRoots(fcomplex a[], int m, fcomplex roots[], int polish, const T eps) {
        void laguer(fcomplex a[], int m, fcomplex *x, int *its);
        int i,its,j,jj;
        fcomplex x,b,c,ad[MAXM];
        for (j=0;j<=m;j++) ad[j]=a[j];
        for (j=m;j>=1;j--) {
            x=std::complex(0.0,0.0);
            if (fabs(x.i) <= 2.0*eps*std::abs(x.r)) x.i=0.0;
            roots[j]=x;
            b=ad[j];
            for (jj=j-1;jj>=0;jj--) {
                c=ad[jj];
                ad[jj]=b;
                b=x*b+c;
            }
        }if (polish)
            for (j=1;j<=m;j++)
                laguer(a,x,eps);
        for (j=2;j<=m;j++) {
            for (i=j-1;i>=1;i--) {
                if (roots[i].r <= x.r) break;
                roots[i+1]=roots[i];
            }
            roots[i+1]=x;
        }
    }
}
#endif