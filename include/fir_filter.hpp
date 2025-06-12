#ifndef FIR_FILTER_H
#define FIR_FILTER_H

#include <array>
#include <algorithm>
#include <memory>
#include <numeric>
#include <numbers>
#include <cmath>
#include <iostream>

class FIR_FILTER{
private:
    size_t _num_taps;
    double* _taps;
    std::allocator<double> _alloc;
public:
    constexpr FIR_FILTER(): _num_taps(0){}

    constexpr size_t getNumTaps() const{return this->_num_taps;}
    constexpr double* getTaps() const{return this->_taps;}

    /**
     * @brief Constructs the taps for a root-raised cosine filter, and automatically adapts the filters variables to be the correct size.
     * @param sps The samples per symbol of the filter.
     * @param beta The roll-off factor.
     * @param span The filter span in symbols.
     */
    constexpr void designTaps(const size_t& sps, const double& beta, const size_t& span){
        if(sps < 1) throw std::logic_error("Samples per symbol cannot be less than 1");
        if(span <= 0) throw std::logic_error("Span must be greater than 0");
        if(beta <= 0 || beta > 1) throw std::logic_error("Beta must be in range (0, 1]");

        double PI = std::numbers::pi_v<double>;
        double eps = 1e-5;
        size_t n = 2*span*sps + 1;
        this->_num_taps = n;
        this->_taps = std::allocator_traits<std::allocator<double>>::allocate(this->_alloc, n);

        for(size_t i = 0; i < n; ++i){
            double nd = static_cast<double>(i);
            double t = nd/static_cast<double>(sps) - static_cast<double>(span);
            if(std::fabs(t) < eps){
                this->_taps[i] = 1.0 + beta*4.0/PI - beta;
            } else if(std::fabs(t) < 1.0/4.0/beta + eps && std::fabs(t) > 1.0/4.0/beta - eps){
                double a1 = (1.0 + 2.0/PI) * std::sin(PI/4.0/beta);
                double a2 = (1.0 - 2.0/PI) * std::sin(PI/4.0/beta);
                this->_taps[i] = beta/std::sqrt(2.0)*(a1 + a2);
            } else{
                double a1 = std::sin(PI*t*(1.0 - beta));
                double a2 = 4.0*beta*t*std::cos(PI*t*(1.0 + beta));
                double a3 = PI*t*(1.0 - 16.0*std::pow(beta, 2.0)*std::pow(t, 2.0));
                this->_taps[i] = (a1 + a2)/a3;
            }
            // std::cout<<_taps[i]<<",";
        }
        // std::cout<<std::endl;
        double tp = std::sqrtl(std::accumulate(this->_taps, this->_taps + n, 0.0));
        for(size_t i = 0; i < n; ++i){
            this->_taps[i] = this->_taps[i]/tp;
        }
        // std::cout<<tp;
    }
};
#endif