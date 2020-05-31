//
//  Gauss_formula.cpp
//  Integral(Task4)
//
//  Created by Артем Белов on 01.06.2020.
//  Copyright © 2020 Артем Белов. All rights reserved.
//
#include "../LUP/LUP.hpp"
#include "../Cardano/Cardano_formula.hpp"

#include <cmath>

#include "Gauss_formula.hpp"

double Gauss_formula() {
    
    const int N = 3;

    // Найдено аналитически
    double *Mu = new double[6];
    Mu[0] = 2.348811462565555; Mu[1] = 2.296615652286315; Mu[2] = 3.248070708233503;
    Mu[3] = 5.265293569136416; Mu[4] = 9.170386299579257; Mu[5] = 16.69642746957879;
    
    int     *p   = new int[N+1];
    double  *b   = new double[N];
    double  *abc = new double[N];
    double  *x   = new double[N];
    double **A   = new double *[N];
    for (int i = 0; i < N; ++i) {
        A[i] = new double[N];
    }
    b[0]    = Mu[3]; b[1]    = Mu[4]; b[2]    = Mu[5];
    
    A[0][0] = Mu[0]; A[0][1] = Mu[1]; A[0][2] = Mu[2];
    A[1][0] = Mu[1]; A[1][1] = Mu[2]; A[1][2] = Mu[3];
    A[2][0] = Mu[2]; A[2][1] = Mu[3]; A[2][2] = Mu[4];
    
    LUPDecompose(A, p, N);
    LUPSolve(A, p, b, abc, N);
    Cardano_formula(x, abc[2], abc[1], abc[0]);
    
    b[0]    = Mu[0];        b[1]    = Mu[1];        b[2]    = Mu[2];
 
    A[0][0] = 1;            A[0][1] = 1;            A[0][2] = 1;
    A[1][0] = x[0];         A[1][1] = x[1];         A[1][2] = x[2];
    A[2][0] = pow(x[0], 2); A[2][1] = pow(x[1], 2); A[2][2] = pow(x[2], 2);

    LUPDecompose(A, p, N);
    LUPSolve(A, p, b, x, N);
    
//    std::cout << x[0] + x[1] + x[2] << std::endl;
    
//    std::cout << std::fixed;
//    std::cout << std::setprecision(54) << x[0] << " " << std::setprecision(50) << x[1] << " " << std::setprecision(49) << x[2] << std::endl;
    
    return x[0] + x[1] + x[2];
}
