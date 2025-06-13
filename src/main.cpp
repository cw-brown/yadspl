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
    b = {1,5,6,10, 30, 12, 19};

    std::cout<<"Before resize: ";
    for(auto&& v:b)
        std::cout<<v<<" ";
    std::cout<<"\n";

    b.resize(3);

    std::cout<<"After resize to degree 3: ";
    for(auto&& v:b)
        std::cout<<v<<" ";
    std::cout<<"\n";

    b.resize(12, 2.5);
    std::cout<<"After resize to degree 12: ";
    for(auto&& v:b)
        std::cout<<v<<" ";
    std::cout<<"\n";


    return 0;
}