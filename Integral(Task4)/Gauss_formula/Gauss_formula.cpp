//
//  Gauss_formula.cpp
//  Integral(Task4)
//
//  Created by Артем Белов on 01.06.2020.
//  Copyright © 2020 Артем Белов. All rights reserved.
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

double Gauss_formula(const double &a, const double &b) {
    
    std::cout << "[ " << a << " ; " << b << " ]\n";
    
    const int N = 3;
    
    double mu[2*N];
    int    *p      = new int[N+1];
    double *b_vect = new double[N];
    double *a_vect = new double[N];
    double *x_vect = new double[N];
    double *A_vect = new double[N];
    double **M     = new double*[N];
    for (int i = 0; i < N; ++i) {
        M[i] = new double[N];
    }
    
    mu[0] = (pow(b - LOWER_LIMIT, 1 - ALPHA_DEGREE) - pow(a - LOWER_LIMIT, 1 - ALPHA_DEGREE))/(1 - ALPHA_DEGREE);
    
    mu[1] = (pow(b - LOWER_LIMIT, 2 - ALPHA_DEGREE) - pow(a - LOWER_LIMIT, 2 - ALPHA_DEGREE))/(2 - ALPHA_DEGREE) + LOWER_LIMIT*mu[0];
    
    mu[2] = (pow(b - LOWER_LIMIT, 3 - ALPHA_DEGREE) - pow(a - LOWER_LIMIT, 3 - ALPHA_DEGREE))/(3 - ALPHA_DEGREE) + 2*LOWER_LIMIT*mu[1]
            - LOWER_LIMIT*LOWER_LIMIT*mu[0];

    mu[3] = (pow(b - LOWER_LIMIT, 4 - ALPHA_DEGREE) - pow(a - LOWER_LIMIT, 4 - ALPHA_DEGREE))/(4 - ALPHA_DEGREE) + 3*LOWER_LIMIT*mu[2]
            -3*LOWER_LIMIT*LOWER_LIMIT*mu[1] + LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*mu[0];
    
    mu[4] = (pow(b - LOWER_LIMIT, 5 - ALPHA_DEGREE) - pow(a - LOWER_LIMIT, 5 - ALPHA_DEGREE))/(5 - ALPHA_DEGREE)
            -LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*mu[0] + 4*LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*mu[1]
            -6*LOWER_LIMIT*LOWER_LIMIT*mu[2] + 4*LOWER_LIMIT*mu[3];
    
    mu[5] = (pow(b - LOWER_LIMIT, 6 - ALPHA_DEGREE) - pow(a - LOWER_LIMIT, 6 - ALPHA_DEGREE))/(6 - ALPHA_DEGREE)
            + LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*mu[0] - 5*LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*mu[1]
            + 10*LOWER_LIMIT*LOWER_LIMIT*LOWER_LIMIT*mu[2] - 10*LOWER_LIMIT*LOWER_LIMIT*mu[3] + 5*LOWER_LIMIT*mu[4];
    
    M[0][0]   = mu[0]; M[0][1]    = mu[1];  M[0][2]   = mu[2];
    M[1][0]   = mu[1]; M[1][1]    = mu[2];  M[1][2]   = mu[3];
    M[2][0]   = mu[2]; M[2][1]    = mu[3];  M[2][2]   = mu[4];
    
    b_vect[0] = -mu[3]; b_vect[1] = -mu[4]; b_vect[2] = -mu[5];
    
    LUPDecompose(M, p, N);
    LUPSolve(M, p, b_vect, a_vect, N);
    
    Cardano_formula(x_vect, a_vect[2], a_vect[1], a_vect[0]);
    
    std::cout << mu[0] << " " << mu[1] << " " << mu[2] << "\n";
    
    std::cout << x_vect[0] << " " << x_vect[1] << " " << x_vect[2] << "\n";
    
    M[0][0]   = 1;                   M[0][1]   = 1;                    M[0][2]   = 1;
    M[0][0]   = x_vect[0];           M[0][1]   = x_vect[1];            M[0][2]   = x_vect[2];
    M[0][0]   = x_vect[0]*x_vect[0]; M[0][1]   = x_vect[1]*x_vect[1];  M[0][2]   = x_vect[2]*x_vect[2];
    
    b_vect[0] = mu[0];               b_vect[1] = mu[1];                b_vect[2] = mu[2];
    
    LUPDecompose(M, p, N);
    LUPSolve(M, p, b_vect, A_vect, N);
    
    for (int i = 0; i < N; ++i) {
        if ((x_vect[i] < a || x_vect[i] > b)&&(x_vect[i] != a && x_vect[i] != b)) {
            std::cout << "!!! " << i << "\n";
        }
    }
    
    double Q, Sn = 0.0;
    for(int i = 0; i < N; ++i) {
        
        Q = (a+b)/2 + (b-a)/2*x_vect[i];
        Sn += A_vect[i]*f_div_p(Q);
    }
    std::cout << (b-a)/2*Sn << "\n";
    return (b-a)/2*Sn;
    
}

double Gauss() {
    
    bool h_opt_flag      = true;
    const double fac     = 0.95;
    
    const double a = LOWER_LIMIT;
    
    const int L = 2;
    int       N = 4;
    unsigned int       m = 0;
    double *S_h = new double[N];
    double    h = 0;
    double  R_h = 0;
    
    do {
        S_h[0] = 0; S_h[1] = 0; S_h[2] = 0;
        h = (HIGHER_LIMIT-LOWER_LIMIT) / N;
        if  (LOWER_LIMIT == h) {
            N++;
            continue;
        }
        
        for (int i = 0; i < N; ++i) {
            
            S_h[0] += Gauss_formula(a + i*h, a + (i+1)*h);
            
        }
        std::cout << "\n";
        for (int i = 0; i < N*L; ++i) {
            S_h[1] += Gauss_formula(a + i*h/L, a + (i+1)*h/L);
        }
        std::cout << "\n";
        for (int i = 0; i < N*L*L; ++i) {
            S_h[2] += Gauss_formula(a + i*h/(L*L), a + (i+1)*h/(L*L));
        }
        m   =  -(log(fabs(S_h[2] - S_h[1])) - log(fabs(S_h[1] - S_h[0])))/(log(L)) + 4;
        

        std::cout << " " << S_h[0] << " " << S_h[1] << " " << S_h[2] << "\n\n";
        
        R_h = fabs(Richardson_extrapolation(S_h, h, m, L));
        
        if (h_opt_flag) {
                    h_opt_flag = false;
                    h = h*pow( MyExp * ((1 - 1/pow(L, m))/(fabs(S_h[1] - S_h[0])))  , 1/m);
                    h *= fac;
                    N = (HIGHER_LIMIT-LOWER_LIMIT) / h;
                    N++;
//                    std::cout << "h_opt = " << h << "\n";
                }
        if (!h_opt_flag) {
             N*=2;
        }
       
    } while (R_h > MyExp);
    
    return S_h[0];
}
