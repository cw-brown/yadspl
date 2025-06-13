#include <iostream>
#include <math.h>
#include <complex>
#include <vector>
#include <typeinfo>
#include <memory>

#include "implot.h"

#include "helper_funcs.h"
#include "fir_filter.hpp"
#include "polynomial.hpp"



int main(){
    // std::vector<double> a = {1,2,3,4};
    polynomial<double> a{1,2,3,4, 20};

    polynomial<double> b;
    b = {1,5,6,10};

    std::cout<<"Front: "<<b.front()<<", Back: "<<b.back()<<"\n";

    std::cout<<"Degree: "<<b.degree()<<"\n";

    std::cout<<"Coefficients: ";
    for(auto&& v:b)
        std::cout<<v<<" ";

    return 0;
}