#ifndef FIR_FILTER_H
#define FIR_FILTER_H

#include <array>
#include <algorithm>
#include <memory>
#include <numeric>
#include <numbers>
#include <cmath>
#include <iostream>

enum class filter_type{
    None,
    RRC,
    Lowpass
};

/**
 * @brief A filter class which holds its own taps. Automatically updates its type depending on which kind of filter was designed.
 * Currently has root-raised cosine, raised cosine, and pole zero filters. All data is double precision.
 * 
 * Approach was based on "Digital Filter Designer's Handbook With C++ Algorithms" by Rorabaugh.
 */
class fir_filter{
public:
    enum class filter_type{
    None,
    RRC,
    RC,
    Lowpass
    };
private:
    size_t _num_taps;
    double* _taps;
    filter_type _type;
    std::allocator<double> _alloc;
public:
    constexpr fir_filter(): _num_taps(0), _taps(nullptr){}

    constexpr size_t getNumTaps() const{return _num_taps;}
    constexpr double* getTaps() const{return _taps;}
    constexpr filter_type getType() const{return _type;}

    /**
     * @brief Designs the filter for a root-raised cosine filter.
     * @param sps The samples per symbol of the filter.
     * @param beta The roll-off factor.
     * @param span The filter span in symbols.
     */
    constexpr void design_rrc(const size_t& sps, const double& beta, const size_t& span){
        if(sps < 1) throw std::logic_error("Samples per symbol cannot be less than 1");
        if(span <= 0) throw std::logic_error("Span must be greater than 0");
        if(beta <= 0 || beta > 1) throw std::logic_error("Beta must be in range (0, 1]");

        _type = filter_type::RRC;
        double PI = std::numbers::pi_v<double>;
        double eps = 1e-5;
        size_t n = 2*span*sps + 1;
        _num_taps = n;
        _taps = std::allocator_traits<std::allocator<double>>::allocate(_alloc, n);

        for(size_t i = 0; i < n; ++i){
            double t = static_cast<double>(i)/static_cast<double>(sps) - static_cast<double>(span);
            if(std::fabs(t) < eps){
                _taps[i] = 1.0 + beta*4.0/PI - beta;
            } else if(std::fabs(t) < 1.0/4.0/beta + eps && std::fabs(t) > 1.0/4.0/beta - eps){
                double a1 = (1.0 + 2.0/PI) * std::sin(PI/4.0/beta);
                double a2 = (1.0 - 2.0/PI) * std::sin(PI/4.0/beta);
                _taps[i] = beta/std::sqrt(2.0)*(a1 + a2);
            } else{
                double a1 = std::sin(PI*t*(1.0 - beta));
                double a2 = 4.0*beta*t*std::cos(PI*t*(1.0 + beta));
                double a3 = PI*t*(1.0 - 16.0*std::pow(beta, 2.0)*std::pow(t, 2.0));
                _taps[i] = (a1 + a2)/a3;
            }
        }
        double tp = std::sqrtl(std::accumulate(_taps, _taps + n, 0.0));
        for(size_t i = 0; i < n; ++i){
            _taps[i] = _taps[i]/tp;
        }
    }

    constexpr void design_rc(const size_t& sps, const double& beta, const size_t& span){
        if(sps < 1) throw std::logic_error("Samples per symbol cannot be less than 1");
        if(span <= 0) throw std::logic_error("Span must be greater than 0");
        if(beta <= 0 || beta > 1) throw std::logic_error("Beta must be in range (0, 1]");

        _type = filter_type::RC;
        double PI = std::numbers::pi_v<double>;
        double eps = 1e-5;
        size_t n = 2*span*sps + 1;
        _num_taps = n;
        _taps = std::allocator_traits<std::allocator<double>>::allocate(_alloc, n);

        for(size_t i = 0; i < n; ++i){
            double t = static_cast<double>(i)/static_cast<double>(sps) - static_cast<double>(span);
            if(std::abs(t) < 1.0/2.0/beta + eps && std::abs(t) > 1.0/2.0/beta - eps){
                _taps[i] = beta/2.0 * std::sin(PI/2.0/beta);
            } else if(std::abs(t) < eps){
                _taps[i] = 1;
            } else{
                double a1 = std::cos(PI*beta*t)/(1 - std::pow(2*beta*t, 2));
                _taps[i] = a1 * std::sin(PI*t)/PI/t;
            }
        }

        double tp = std::sqrtl(std::accumulate(_taps, _taps + n, 0.0));
        for(size_t i = 0; i < n; ++i){
            _taps[i] = _taps[i]/tp;
        }

    }
};
#endif