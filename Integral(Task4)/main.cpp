//
//  main.cpp
//  Integral(Task4)
//
//  Created by Артем Белов on 28.05.2020.
//  Copyright © 2020 Артем Белов. All rights reserved.
//

// Ответ Wolfram: 3.57886
// Совпадает с точным решением, если разбить на 6682 узлов, при этом h = = 0.000329292




#include <iostream>
#include <cmath>

#include "LUP/LUP.hpp"
#include "constants.cpp"


#include "Gauss_formula/Gauss_formula.hpp"
#include "Newton_Cotes/Newton_Cotes.hpp"



int main() {

    double Newton_Cotes_result = Newton_Cotes();
    double Gauss_result = Gauss_formula();
    
    std::cout << "Результат с помощью формул Ньютона-Котеса: " << Newton_Cotes_result << "\n";
    std::cout << "Результат с помощью формулы Гаусса: " << Gauss_result << "\n";
    
    
    
    return 0;
}
