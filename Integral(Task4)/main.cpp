//
//  main.cpp
//  Integral(Task4)
//
//  Created by Артем Белов on 28.05.2020.
//  Copyright © 2020 Артем Белов. All rights reserved.
//

/*
https://www.wolframalpha.com/input/?i=integrate+%28%282.5*cos%282*x%29*exp%282*x%2F3%29%2B4*sin%283.5*x%29*exp%28-3*x%29%2B3*x%29%2F%28x+-+0.1%29%5E%280.2%29%29+dx+from+x%3D0.1+to+2.3
 */
// Ответ Wolfram: 3.57886
// Совпадает с точным решением, если разбить на 6682 узлов, при этом h = 0.000329292

#include <iostream>

#include "Gauss_formula/Gauss_formula.hpp"
#include "Newton_Cotes/Newton_Cotes.hpp"


int main() {

    double Newton_Cotes_result = Newton_Cotes();
    double Gauss_result = Gauss();
    
    std::cout << "Результат с помощью формул Ньютона-Котеса: " << Newton_Cotes_result << "\n";
    std::cout << "Результат с помощью формулы Гаусса: " << Gauss_result << "\n";
    
    
    
    return 0;
}
