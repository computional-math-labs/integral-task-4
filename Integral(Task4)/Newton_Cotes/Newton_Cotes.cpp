//
//  Newton_Cotes.cpp
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

#include "Newton_Cotes.hpp"

double Newton_Cotes_formula(const double &h, const int &N) {
    
    double Sn = 0.0;
    
    const unsigned short int M = 3;
    
    double* Z = new double[N+1];
    double temp = LOWER_LIMIT;
    for (int i = 0; i < N+1; ++i) {
        Z[i] = temp;
        temp+=h;
    }

    
    double** Mu = new double*[N+1];
    long double** A = new long double*[N+1];
    for (int i = 0; i < N; ++i) {
        Mu[i] = new double[M];
        A[i] = new long double[M];
    }

    for (int i = 1; i < N+1; ++i) {
        Mu[i-1][0] = (pow(Z[i] - LOWER_LIMIT, 1 - ALPHA_DEGREE) - pow(Z[i-1] - LOWER_LIMIT, 1 - ALPHA_DEGREE))/(1 - ALPHA_DEGREE);
        Mu[i-1][1] = (pow(Z[i] - LOWER_LIMIT, 2 - ALPHA_DEGREE) - pow(Z[i-1] - LOWER_LIMIT, 2 - ALPHA_DEGREE))/(2 - ALPHA_DEGREE)
                      + LOWER_LIMIT*Mu[i-1][0];
        Mu[i-1][2] = (pow(Z[i] - LOWER_LIMIT, 3 - ALPHA_DEGREE) - pow(Z[i-1] - LOWER_LIMIT, 3 - ALPHA_DEGREE))/(3 - ALPHA_DEGREE)
                      + 2*LOWER_LIMIT*Mu[i-1][1] - LOWER_LIMIT*LOWER_LIMIT*Mu[i-1][0];
        
        A[i-1][0]  =  (Mu[i-1][2] - Mu[i-1][1]*(Z[i]/2 + Z[i]) + Mu[i-1][0]*Z[i]/2*Z[i])/((Z[i]/2-Z[i-1])*(Z[i] - Z[i-1]));
        A[i-1][1]  = -(Mu[i-1][2] -  Mu[i-1][1]*(Z[i-1] + Z[i]) + Mu[i-1][0]*Z[i-1]*Z[i])/((Z[i]/2 - Z[i-1])*(Z[i] - Z[i]/2));
        A[i-1][2]  =  (Mu[i-1][2] -  Mu[i-1][1]*(Z[i]/2 + Z[i-1]) + Mu[i-1][0]*Z[i]/2*Z[i-1]  )/((Z[i] - Z[i]/2)*(Z[i] - Z[i-1]));
    }
    
    Sn += A[0][0]*f(Z[0]) + A[N-1][2]*f(Z[N-1]);
    
    for (int i = 1; i < N+1; ++i) {
        Sn += A[i-1][1]*f(Z[i]/2);
    }
    
    for (int i = 1; i < N; ++i) {
        Sn += (A[i-1][2] + A[i][0])*f(Z[i]);
    }
    
    if (Z[1]/2-Z[0] == 0) {
        std::cout << "Z_0 " << Z[1]/2-Z[0] << " h =" << h << " N = " << N << std::endl;
    }
    
    return Sn;
}

double Newton_Cotes() {
    
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
                
        S_h[0] =  Newton_Cotes_formula( h, N);
        S_h[1] =  Newton_Cotes_formula( h/L, L*N);
        S_h[2] =  Newton_Cotes_formula( h/(L*L), L*L*N);
        
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
