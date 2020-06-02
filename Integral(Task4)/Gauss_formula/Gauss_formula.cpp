//
//  Gauss_formula.cpp
//  Integral(Task4)
//
//  Created by Артем Белов on 01.06.2020.
//  Copyright © 2020 Артем Белов. All rights reserved.
//
//
#include "../constants.cpp"
#include "../f/f.hpp"
#include "../LUP/LUP.hpp"
#include "../LUP/Gauss_SLAE.hpp"
#include "../Cardano/Cardano_formula.hpp"
#include "../Richardson_extrapolation/Richardson_extrapolation.hpp"

#include <cmath>
#include <iostream>

#include "Gauss_formula.hpp"

double Gauss_formula(const double &h, const int &N) {
     
    double Sh  = 0;
    double *Z  = new double[N+1];
    int    *p  = new int[3];
    double* x = new double[3];
    
    double* mu = new double[2*N-1];
    double* A = new double[3];
    double* a = new double[3];
    double* c = new double[3];
    double** M = new double* [3];
    for (int i = 0; i < 3; ++i) {
        M[i] = new double[3];
    }

    for (int i = 0; i < N + 1; ++i) {
        Z[i] = LOWER_LIMIT + i * h;
    }


    for (int i = 1; i < N + 1; ++i) {

        mu[0] = (pow(Z[i] - LOWER_LIMIT, 1 - ALPHA_DEGREE) - pow(Z[i-1] - LOWER_LIMIT, 1 - ALPHA_DEGREE))/(1 - ALPHA_DEGREE);
        
        mu[1] = (pow(Z[i] - LOWER_LIMIT, 2 - ALPHA_DEGREE) - pow(Z[i-1] - LOWER_LIMIT, 2 - ALPHA_DEGREE))/(2 - ALPHA_DEGREE) + LOWER_LIMIT*mu[0];
        
        mu[2] = (pow(Z[i] - LOWER_LIMIT, 3 - ALPHA_DEGREE) - pow(Z[i-1] - LOWER_LIMIT, 3 - ALPHA_DEGREE))/(3 - ALPHA_DEGREE) + 2*LOWER_LIMIT*mu[1]
                - LOWER_LIMIT*LOWER_LIMIT*mu[0];

        mu[3] = (pow(Z[i] - LOWER_LIMIT, 4 - ALPHA_DEGREE) - pow(Z[i-1] - LOWER_LIMIT, 4 - ALPHA_DEGREE))/(4 - ALPHA_DEGREE) + 3*LOWER_LIMIT*mu[2]
                -3*LOWER_LIMIT*LOWER_LIMIT*mu[1] + LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*mu[0];
        
        mu[4] = (pow(Z[i] - LOWER_LIMIT, 5 - ALPHA_DEGREE) - pow(Z[i-1] - LOWER_LIMIT, 5 - ALPHA_DEGREE))/(5 - ALPHA_DEGREE)
                -LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*mu[0] + 4*LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*mu[1]
                -6*LOWER_LIMIT*LOWER_LIMIT*mu[2] + 4*LOWER_LIMIT*mu[3];
        
        mu[5] = (pow(Z[i] - LOWER_LIMIT, 6 - ALPHA_DEGREE) - pow(Z[i-1]- LOWER_LIMIT, 6 - ALPHA_DEGREE))/(6 - ALPHA_DEGREE)
                + LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*mu[0] - 5*LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*mu[1]
                + 10*LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*mu[2] - 10*LOWER_LIMIT*LOWER_LIMIT*mu[3] + 5*LOWER_LIMIT*mu[4];

        M[0][0] = mu[0]; M[0][1] = mu[1]; M[0][2] = mu[2];
        M[1][0] = mu[1]; M[1][1] = mu[2]; M[1][2] = mu[3];
        M[2][0] = mu[2]; M[2][1] = mu[3]; M[2][2] = mu[4];

        c[0]    = -mu[3]; c[1]   = -mu[4]; c[2]   = -mu[5];


        LUPDecompose(M, p, 3);
        LUPSolve(M, p, c, a, 3);

        Cardano_formula(x, a[2], a[1], a[0]);


        M[0][0] = 1;         M[0][1] = 1;         M[0][2] = 1;
        M[1][0] = x[0];      M[1][1] = x[1];      M[1][2] = x[2];
        M[2][0] = x[0]*x[0]; M[2][1] = x[1]*x[1]; M[2][2] = x[2] * x[2];

        c[0]    = mu[0];     c[1]    = mu[1];     c[2]    = mu[2];


        LUPDecompose(M, p, 3);
        LUPSolve(M, p, c, A, 3);
        Sh += A[0] * f(x[0]) + A[1] * f(x[1]) + A[2] * f(x[2]);

    }

    return Sh;
}

double Gauss() {
    
    bool h_opt_flag      = true;
        
        const int    L       = 2;
        const double fac     = 0.95;
        
        unsigned int N       = 8;
        double       *S_h    = new double[3];
        double       h       = 0;
        unsigned int m       = 0;
        double       R_h     = 0;
        
        do {
            
            h = (HIGHER_LIMIT-LOWER_LIMIT) / N;
            if (LOWER_LIMIT == h) {
                N++;
                continue;
            }
                    
            S_h[0] =  Gauss_formula(h, N);
            S_h[1] =  Gauss_formula( h/L, L*N);
            S_h[2] =  Gauss_formula( h/(L*L), L*L*N);
            
            m      =  -(log(fabs(S_h[2] - S_h[1])) - log(fabs(S_h[1] - S_h[0])))/(log(L)) + 1;
            
            R_h    =  fabs(Richardson_extrapolation(S_h, h, m, L));
            
            if (h_opt_flag) {
                h_opt_flag = false;
                h = h*pow( MyExp * ((1 - 1/pow(L, m))/(fabs(S_h[1] - S_h[0])))  , 1/m);
                h *= fac;
                N = (HIGHER_LIMIT-LOWER_LIMIT) / h;
                N++;
    //            std::cout << "h_opt = " << h << "\n";
            }
            
            N*=2;
        
        } while (R_h > MyExp);
        
        return S_h[0];
}
