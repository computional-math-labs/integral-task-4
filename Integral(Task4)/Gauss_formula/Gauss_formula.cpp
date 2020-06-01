//
//  Gauss_formula.cpp
//  Integral(Task4)
//
//  Created by Артем Белов on 01.06.2020.
//  Copyright © 2020 Артем Белов. All rights reserved.
//
#include "../f/f.hpp"
#include "../LUP/LUP.hpp"
#include "../Cardano/Cardano_formula.hpp"

#include <cmath>
#include <iostream>

#include "Gauss_formula.hpp"

double Gauss_formula() {
    
    const int N = 3;

    // Найдено аналитически и проверено в Wolfram. Ссылка на калькулятор:
    /*
     https://www.wolframalpha.com/input/?i=integrate+x%5E4*%28x-0.1%29%5E%28-0.2%29+dx+from+x%3D0.1+to+2.3
     */
    double *Mu = new double[6];
    Mu[0] = 2.34881146256555;  Mu[1] = 2.53149679854287;  Mu[2]  = 3.730881953316422;
    Mu[3] = 6.310962062637622; Mu[4] = 11.48080931348324; Mu[5]  = 21.84180247930516;
    
    int     *p   = new int[N+1];
    double  *b   = new double[N];
    double  *a   = new double[N];
    double  *A   = new double[N];
    double  *t   = new double[N];
    double **M   = new double *[N];
    for (int i = 0; i < N; ++i) {
        M[i] = new double[N];
    }
    b[0]    = -Mu[3]; b[1]   = -Mu[4]; b[2]   = -Mu[5];
    
    M[0][0] = Mu[0]; M[0][1] = Mu[1]; M[0][2] = Mu[2];
    M[1][0] = Mu[1]; M[1][1] = Mu[2]; M[1][2] = Mu[3];
    M[2][0] = Mu[2]; M[2][1] = Mu[3]; M[2][2] = Mu[4];
    
    LUPDecompose(M, p, N);
    LUPSolve(M, p, b, a, N);
    
    std::cout << a[2] << " " << a[1] << " " << a[0] << "\n";
    
    Cardano_formula(t, a[2], a[1], a[0]);
    
    std::cout << t[0] << " " << t[1] << " " << t[2];
    
    b[0]    = Mu[0];        b[1]    = Mu[1];        b[2]    = Mu[2];
    
    
 
    M[0][0] = 1;            M[0][1] = 1;            M[0][2] = 1;
    M[1][0] = t[0];         M[1][1] = t[1];         M[1][2] = t[2];
    M[2][0] = pow(t[0], 2); M[2][1] = pow(t[1], 2); M[2][2] = pow(t[2], 2);

    LUPDecompose(M, p, N);
    LUPSolve(M, p, b, A, N);
    
    return A[0]*f(t[0]) + A[1]*f(t[1]) + A[2]*f(t[2]);
}
