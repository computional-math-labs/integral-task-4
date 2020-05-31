//
//  main.cpp
//  Integral(Task4)
//
//  Created by Артем Белов on 28.05.2020.
//  Copyright © 2020 Артем Белов. All rights reserved.
//

// Ответ Wolfram: 3.57886
// Совпадает с точным решением на 20 шаге


#define  MyExp 0.00000001
#define LOWER_LIMIT 0.1
#define HIGHER_LIMIT 2.3
#define ALPHA_DEGREE 0.2

#include <iostream>
#include <cmath>

#include "LUP/LUP.hpp"

void newline() {
    std::cout << std::endl;
}

double f(double x) {
    return (2.5*cos(2*x)*exp(2*x/3)+4*sin(3.5*x)*exp(-3*x)+3*x);
}

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
    
    
//    for (int i = 0; i < N; ++i) {
//        for (int j = 0; j < M; ++j) {
//            std::cout << Mu[i][j] << " ";
//        }
//        newline();
//    }
//    newline();
//
//    for (int i = 0; i < N; ++i) {
//        for (int j = 0; j < M; ++j) {
//            std::cout << A[i][j] << " ";
//        }
//        newline();
//    }
//    newline();
    
    return Sn;
}

double Richardson_extrapolation(double *S_h, const double &h, const int &m, const int &L) {
    
    double R_h = 0;
    
    
    double *Cm = new double[3];
    
    int *p = new int[3];
    double **H = new double*[3];
    for (int i = 0; i < 3; ++i) {
        H[i] = new double[3];
    }
    H[0][0] = 1; H[0][1] = -pow(h, m);       H[0][2] = -pow(h, m+1);
    H[1][0] = 1; H[1][1] = -pow(h/L, m);     H[1][2] = -pow(h/L, m+1);
    H[2][0] = 1; H[2][1] = -pow(h/(L*L), m); H[2][2] = -pow(h/(L*L), m+1);
    
    LUPDecompose(H, p, 3);
    LUPSolve(H, p, S_h, Cm, 3);
    
    R_h = Cm[0] - S_h[0];
    
    std::cout << "J - S_h = R_h = " << R_h  << std::endl;
    return R_h;
}

double integral() {
    
    const int L = 2;
    
    unsigned int N      = 3;
    double *S_h = new double[3];
    double       h      = 0;
    unsigned int m      = 0;
//    double       R_h    = 0;
    
    do {
        h = (HIGHER_LIMIT-LOWER_LIMIT) / N;
        
        if (LOWER_LIMIT == h) {
            N++;
            continue;
        }
                
                S_h[0] = Newton_Cotes_formula( h, N);
                S_h[1]= Newton_Cotes_formula( h/L, L*N);
                S_h[2] = Newton_Cotes_formula( h/(L*L), L*L*N);
                
                m = -(log(fabs(S_h[2] - S_h[1])) - log(fabs(S_h[1] - S_h[0])))/(log(L)) + 1;
        
        Richardson_extrapolation(S_h, h, m, L);
        N*=8;
        std::cout << "N = " << N << " S_h = " << S_h[0] << " m = " << m << "\n\n";
    
    } while (Richardson_extrapolation(S_h, h, m, L) > MyExp);
    
    
    return S_h[0];
}

int main() {

    double Newton_Cotes_result = integral();
    std::cout << Newton_Cotes_result;
    newline();

    return 0;
}
