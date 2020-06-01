//
//  Richardson_extrapolation.cpp
//  Integral(Task4)
//
//  Created by Артем Белов on 01.06.2020.
//  Copyright © 2020 Артем Белов. All rights reserved.
//

#include <cmath>

#include "../LUP/LUP.hpp"
#include "Richardson_extrapolation.hpp"

double Richardson_extrapolation(double *S_h, const double &h, const int &m, const int &L) {
    
    double R_h = 0;
    
    double *Cm = new double[3];
    int     *p = new int[3];
    double **H = new double*[3];
    for (int i = 0; i < 3; ++i) {
        H[i] = new double[3];
    }
    
    H[0][0] = 1; H[0][1] = -pow(h, m);       H[0][2] = -pow(h, m+1);
    H[1][0] = 1; H[1][1] = -pow(h/L, m);     H[1][2] = -pow(h/L, m+1);
    H[2][0] = 1; H[2][1] = -pow(h/(L*L), m); H[2][2] = -pow(h/(L*L), m+1);
    
    if (m == -1) {
        H[0][1] = -1/pow(h, -m);
        H[1][1] = -1/pow(h/L, -m);
        H[2][1] = -1/pow(h/(L*L), -m);
    } else if (m < -1) {
        H[0][2] = -1/pow(h, -(m+1));
        H[1][2] = -1/pow(h/L, -(m+1));
        H[2][2] = -1/pow(h/(L*L), -(m+1));
    } 
    
    LUPDecompose(H, p, 3);
    LUPSolve(H, p, S_h, Cm, 3);
    
    R_h = Cm[0] - S_h[0];

    return R_h;
}
