//
//  Cardano_formula.cpp
//  Integral(Task4)
//
//  Created by Артем Белов on 01.06.2020.
//  Copyright © 2020 Артем Белов. All rights reserved.
//
#include "../constants.cpp"

#include <cmath>

#include "Cardano_formula.hpp"

int Cardano_formula(double *x,double a,double b,double c) {
  double q,r,r2,q3;
  q=(a*a-3.*b)/9.; r=(a*(2.*a*a-9.*b)+27.*c)/54.;
  r2=r*r; q3=q*q*q;
  if(r2<q3) {
    double t=acos(r/sqrt(q3));
    a/=3.; q=-2.*sqrt(q);
    x[0]=q*cos(t/3.)-a;
    x[1]=q*cos((t+M_2PI)/3.)-a;
    x[2]=q*cos((t-M_2PI)/3.)-a;
    return(3);
  }
  else {
    double aa,bb;
    if(r<=0.) r=-r;
    aa=-pow(r+sqrt(r2-q3),1./3.);
    if(aa!=0.) bb=q/aa;
    else bb=0.;
    a/=3.; q=aa+bb; r=aa-bb;
    x[0]=q-a;
    x[1]=(-0.5)*q-a;
    x[2]=(sqrt(3.)*0.5)*fabs(r);
    if(x[2]==0.) return(2);
    return(1);
  }
}
