#include <iostream>
#include <math.h>
#include <complex>
#include <vector>
#include <typeinfo>
#include <memory>

// #include "implot.h"

// #include "helper_funcs.h"
// #include "fir_filter.hpp"
#include "polynomial.hpp"

int main(){
    poly<double, 5> a = {1,2,3,4,5,2};
    poly<double, 5> b = {2,4,5,6,2,1};

    std::cout<<std::boolalpha<<(a==b);

    // poly<double, 20> a;

    // std::cout<<"Degree: "<<a.degree()<<"\n";
    // std::cout<<"Coefficients: ";
    // for(auto&& v:a)
    //     std::cout<<v<<" ";
    // std::cout<<"\n";


    return 0;
}