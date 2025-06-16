#ifndef FIR_FILTER_H
#define FIR_FILTER_H

#include <array>
#include <algorithm>
#include <memory>
#include <numeric>
#include <numbers>
#include <cmath>
#include <iostream>
#include <complex>

enum class filter_type{
    None,
    Root_Raised_Cosine
};

/**
 * @brief A filter class which holds its own taps. Automatically updates its type depending on which kind of filter was designed.
 * Currently has root-raised cosine, raised cosine, and pole zero filters. All data is double precision.
 * 
 * Approach was based on "Digital Filter Designer's Handbook With C++ Algorithms" by Rorabaugh.
 */
template<class __T, filter_type __F>
requires std::is_floating_point_v<__T>
class fir_filter{
public:
    typedef __T value_type;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef __T* pointer;
    typedef __T& reference;
    typedef std::allocator_traits<std::allocator<__T>>::allocator_type allocator_type;
private:
    size_type _n;
    pointer _taps;
    pointer _resp;
    pointer _freqs;
    std::allocator<__T> _alloc;
    value_type PI = std::numbers::pi_v<__T>;
public:
    constexpr fir_filter(): _n(0), _taps(nullptr){}
    constexpr fir_filter(const fir_filter& other): _n(other._n), _taps(other._taps){}
    constexpr fir_filter(fir_filter&& other) noexcept(true): _n(other._n), _taps(std::move(other._taps)){other._taps = nullptr;}

    constexpr size_type getNumTaps() const noexcept{return _n;}
    constexpr pointer getTaps() const noexcept{return _taps;}
    constexpr pointer getMagnitudeResponse() const noexcept{return _resp;}
    constexpr pointer getFrequencyPoints() const noexcept{return _freqs;}

    constexpr void design(const size_type& sps, const value_type& beta, const size_type& span) 
    requires(__F == filter_type::Root_Raised_Cosine){
        if(sps < 1.0) throw std::logic_error("Samples per symbol cannot be less than 1");
        if(span <= 0.0) throw std::logic_error("Span must be greater than 0 symbols");
        if(beta <= 0.0 || beta > 1.0) throw std::logic_error("Beta must be in range (0, 1]");

        value_type eps = static_cast<__T>(1e-5);
        value_type power = 0.0;
        size_type n = 2*span*sps + 1;
        _n = n;
        _taps = _alloc.allocate(n);
        for(size_type i = 0; i < n; ++i){
            value_type t = static_cast<value_type>(i)/static_cast<value_type>(sps) - static_cast<value_type>(span);
            if(std::abs(t) < eps){
                _taps[i] = 1.0 + beta*4.0/PI - beta;
            } else if(std::abs(t) < 1.0/4.0/beta + eps && std::abs(t) > 1.0/4.0/beta - eps){
                value_type a1 = (1.0 + 2.0/PI) * std::sin(PI/4.0/beta);
                value_type a2 = (1.0 - 2.0/PI) * std::sin(PI/4.0/beta);
                _taps[i] = beta/std::sqrt(2.0)*(a1 + a2);
            } else{
                value_type a1 = std::sin(PI*t*(1.0 - beta));
                value_type a2 = 4.0*beta*t*std::cos(PI*t*(1.0 + beta));
                value_type a3 = PI*t*(1.0 - 16.0*std::pow(beta, 2.0)*std::pow(t, 2.0));
                _taps[i] = (a1 + a2)/a3;
            }
            power += _taps[i];
        }
        power = std::sqrt(power);
        for(size_type i = 0; i < n; ++i){
            _taps[i] /= power;
        }
    }

    constexpr void response(const size_type& points){
        _resp = _alloc.allocate(points);
        _freqs = _alloc.allocate(points);
        for(size_type idx = 0; idx < points; ++idx){
            __T lambda = idx*PI/points;
            std::complex<__T> work(0.0, 0.0);
            for(size_type tap_idx = 0; tap_idx < _n; ++tap_idx){
                work += _taps[tap_idx] * std::complex<__T>(std::cos(tap_idx*lambda), -1.0*std::sin(tap_idx*lambda));
            }
            _freqs[idx] = lambda;
            _resp[idx] = 20.0 * std::log10(std::pow(std::sqrt(std::abs(work)), 2.0));
        }
    }
};





// /**
//  * @brief A filter class which holds its own taps. Automatically updates its type depending on which kind of filter was designed.
//  * Currently has root-raised cosine, raised cosine, and pole zero filters. All data is double precision.
//  * 
//  * Approach was based on "Digital Filter Designer's Handbook With C++ Algorithms" by Rorabaugh.
//  */
// class fir_filter{
// public:
// private:
//     size_t _num_taps;
//     double* _taps;
//     filter_type _type;
//     std::allocator<double> _alloc;
// public:
//     constexpr fir_filter(): _num_taps(0), _taps(nullptr){}

//     constexpr size_t getNumTaps() const{return _num_taps;}
//     constexpr double* getTaps() const{return _taps;}
//     constexpr filter_type getType() const{return _type;}

//     /**
//      * @brief Designs the filter for a root-raised cosine filter.
//      * @param sps The samples per symbol of the filter.
//      * @param beta The roll-off factor.
//      * @param span The filter span in symbols.
//      */
//     constexpr void design_rrc(const size_t& sps, const double& beta, const size_t& span){
//         if(sps < 1) throw std::logic_error("Samples per symbol cannot be less than 1");
//         if(span <= 0) throw std::logic_error("Span must be greater than 0");
//         if(beta <= 0 || beta > 1) throw std::logic_error("Beta must be in range (0, 1]");

//         _type = filter_type::RRC;
//         double PI = std::numbers::pi_v<double>;
//         double eps = 1e-5;
//         size_t n = 2*span*sps + 1;
//         _num_taps = n;
//         _taps = std::allocator_traits<std::allocator<double>>::allocate(_alloc, n);

//         for(size_t i = 0; i < n; ++i){
//             double t = static_cast<double>(i)/static_cast<double>(sps) - static_cast<double>(span);
//             if(std::fabs(t) < eps){
//                 _taps[i] = 1.0 + beta*4.0/PI - beta;
//             } else if(std::fabs(t) < 1.0/4.0/beta + eps && std::fabs(t) > 1.0/4.0/beta - eps){
//                 double a1 = (1.0 + 2.0/PI) * std::sin(PI/4.0/beta);
//                 double a2 = (1.0 - 2.0/PI) * std::sin(PI/4.0/beta);
//                 _taps[i] = beta/std::sqrt(2.0)*(a1 + a2);
//             } else{
//                 double a1 = std::sin(PI*t*(1.0 - beta));
//                 double a2 = 4.0*beta*t*std::cos(PI*t*(1.0 + beta));
//                 double a3 = PI*t*(1.0 - 16.0*std::pow(beta, 2.0)*std::pow(t, 2.0));
//                 _taps[i] = (a1 + a2)/a3;
//             }
//         }
//         double tp = std::sqrtl(std::accumulate(_taps, _taps + n, 0.0));
//         for(size_t i = 0; i < n; ++i){
//             _taps[i] = _taps[i]/tp;
//         }
//     }

//     constexpr void design_rc(const size_t& sps, const double& beta, const size_t& span){
//         if(sps < 1) throw std::logic_error("Samples per symbol cannot be less than 1");
//         if(span <= 0) throw std::logic_error("Span must be greater than 0");
//         if(beta <= 0 || beta > 1) throw std::logic_error("Beta must be in range (0, 1]");

//         _type = filter_type::RC;
//         double PI = std::numbers::pi_v<double>;
//         double eps = 1e-5;
//         size_t n = 2*span*sps + 1;
//         _num_taps = n;
//         _taps = std::allocator_traits<std::allocator<double>>::allocate(_alloc, n);

//         for(size_t i = 0; i < n; ++i){
//             double t = static_cast<double>(i)/static_cast<double>(sps) - static_cast<double>(span);
//             if(std::abs(t) < 1.0/2.0/beta + eps && std::abs(t) > 1.0/2.0/beta - eps){
//                 _taps[i] = beta/2.0 * std::sin(PI/2.0/beta);
//             } else if(std::abs(t) < eps){
//                 _taps[i] = 1;
//             } else{
//                 double a1 = std::cos(PI*beta*t)/(1 - std::pow(2*beta*t, 2));
//                 _taps[i] = a1 * std::sin(PI*t)/PI/t;
//             }
//         }

//         double tp = std::sqrtl(std::accumulate(_taps, _taps + n, 0.0));
//         for(size_t i = 0; i < n; ++i){
//             _taps[i] = _taps[i]/tp;
//         }

//     }

//     constexpr void butterworth(const size_t){

//     }
// };
#endif