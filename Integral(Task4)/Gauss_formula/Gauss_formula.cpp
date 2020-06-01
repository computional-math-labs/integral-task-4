//
//  Gauss_formula.cpp
//  Integral(Task4)
//
//  Created by Артем Белов on 01.06.2020.
//  Copyright © 2020 Артем Белов. All rights reserved.
//
#include "../constants.cpp"
#include "../f/f.hpp"
#include "../Richardson_extrapolation/Richardson_extrapolation.hpp"

#include <cmath>
#include <iostream>

#include "Gauss_formula.hpp"

double Gauss_formula(const double &a, const double &b) {
    
    const int N = 3;
    
    const double t[N]={-sqrt(0.6), 0, sqrt(0.6)};
    const double A[N]={5.0/9.0, 8.0/9.0, 5.0/9.0};

    
    double Q, Sn = 0.0;
    for(int i = 0; i < N; ++i) {
        
        Q = (a+b)/2 + (b-a)/2*t[i];
        Sn += A[i]*f_div_p(Q);
    }
    return (b-a)/2*Sn;
    
}

double Gauss() {
    
    bool h_opt_flag      = true;
    const double fac     = 0.95;
    
    const double a = LOWER_LIMIT;
    
    const int L = 2;
    int       N = 2;
    int       m = 0;
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
        
        for (int i = 0; i < N*L; ++i) {
            S_h[1] += Gauss_formula(a + i*h/L, a + (i+1)*h/L);
        }
        
        for (int i = 0; i < N*L*L; ++i) {
            S_h[2] += Gauss_formula(a + i*h/(L*L), a + (i+1)*h/(L*L));
        }

        m   =  -(log(fabs(S_h[2] - S_h[1])) - log(fabs(S_h[1] - S_h[0])))/(log(L)) + 4;
        

//        std::cout << " " << S_h[0] << " " << S_h[1] << " " << S_h[2] << "\n";
        
        R_h = fabs(Richardson_extrapolation(S_h, h, m, L));
        
        if (h_opt_flag) {
                    h_opt_flag = false;
                    h = h*pow( MyExp * ((1 - 1/pow(L, m))/(fabs(S_h[1] - S_h[0])))  , 1/m);
                    h *= fac;
                    N = (HIGHER_LIMIT-LOWER_LIMIT) / h;
                    N++;
//                    std::cout << "h_opt = " << h << "\n";
                }
        
        N*=2;
        
    } while (R_h > MyExp);
    
    return S_h[0];
}
