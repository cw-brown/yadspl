/**
 * @file yadpsl_math.hpp
 * @author cw-brown (https://github.com/cw-brown)
 * @brief Implementation of various math functions.
 * @version 0.1
 * @date 2025-06-18
 */
#ifndef YADSPL_MATH_H
#define YADSPL_MATH_H

#include <complex>
#include <utility>
#include <numeric>
#include <cmath>

namespace yadspl{

    /**
     * @brief The most extremely slow version of the DFT. Do not use this.
     */
    template<typename _Tp> 
    requires std::is_floating_point_v<_Tp>
    std::tuple<_Tp*, _Tp*, _Tp*> DO_NOT_USE_slow_dft(const _Tp* x, const std::size_t& n){
        using namespace std::literals::complex_literals;
        const _Tp PI = std::numbers::pi_v<_Tp>;
        const _Tp d_theta = 2.0 * PI / static_cast<_Tp>(n);
        std::tuple<_Tp*, _Tp*, _Tp*> out;
        std::get<0>(out) = new _Tp[n];
        std::get<1>(out) = new _Tp[n];
        std::get<2>(out) = new _Tp[n];
        for(std::size_t i = 0; i < n; ++i){
            _Tp theta = i*d_theta;
            std::complex<_Tp> work(_Tp{}, _Tp{});
            for(std::size_t j = 0; j < n; ++j){
                work += x[j]*std::exp(static_cast<_Tp>(j)*-1.0i*theta);
            }
            std::get<0>(out)[i] = theta;
            std::get<1>(out)[i] = std::abs(work);
            std::get<2>(out)[i] = std::arg(work);
        }
        return out;
    }



}
#endif