#ifndef FIR_FILTER_H
#define FIR_FILTER_H

#include <algorithm>
#include <memory>
#include <numeric>
#include <numbers>
#include <cmath>
#include <iostream>
#include <complex>
#include <utility>

#include "basic_iterators.hpp"

enum class filter_type{
    /* Basic Root Raised Cosine */
    Root_Raised_Cosine,
    /* Ideal lowpass filter based on the Fourier series method */
    Lowpass_Ideal,
    /* Ideal highpass filter based on the Fourier series method */
    Highpass_Ideal,
    /* Ideal bandpass filter based on the Fourier series method */
    Bandpass_Ideal,
    /* Ideal bandstop filter based on the Fourier series method */
    Bandstop_Ideal,
    /* Filter designed using frequency samples */
    Frequency_Sampled,
    /* Filter designed using frequency samples and optimized for stopband ripple */
    Frequency_Sampled_Optimized
};

enum class window_type{
    None,
    Rectangular,
    Triangular,
    Von_Hann,
    Hamming,
    Blackman_Harris,
    Gaussian,
    Tukey
};

/**
 * @brief A filter class which holds its own taps. Automatically updates its type depending on which kind of filter was designed.
 * All FIR filters must have a type, and different types are not interoperable. Filters can also be provided with a window that will modify the taps. 
 * 
 * HOW TO USE:
 * 
 * 1) Call design() and provide the arguments necessary to construct the specific filter type.
 * 
 * 2) Call synthesize(). If your filter does not have a window, this isn't strictly necessary, but doing it for a None window will optimize to a no-op, and it will be easy to implement a window without changing any other code.
 * 
 * Approach was based on "Digital Filter Designer's Handbook With C++ Algorithms" by Rorabaugh.
 */
template<class __T, filter_type __F, window_type __W = window_type::None>
requires std::is_floating_point_v<__T>
class fir_filter{
public:
    typedef __T value_type;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef __T* pointer;
    typedef __T& reference;
    typedef __iterator<__T> iterator;
    typedef __const_iterator<__T> const_iterator;
    typedef __reverse_iterator<__T> reverse_iterator;
    typedef __const_reverse_iterator<__T> const_reverse_iterator;
    typedef std::allocator_traits<std::allocator<__T>>::allocator_type allocator_type;
private:
    size_type _n;
    pointer _taps;
    std::allocator<__T> _alloc;
    static constexpr value_type PI = std::numbers::pi_v<__T>;
    static constexpr value_type eps = static_cast<__T>(1e-5);
public:
    constexpr fir_filter(): _n(0), _taps(nullptr){}
    constexpr fir_filter(const fir_filter& other): _n(other._n), _taps(_alloc.allocate(_n)){
        std::uninitialized_copy_n(other._taps, other._n, _taps);
    }
    constexpr fir_filter(fir_filter&& other) noexcept(true): _n(other._n), _taps(std::move(other._taps)){other._taps = nullptr;}
    constexpr fir_filter& operator=(const fir_filter& other){
        if(this != &other){
            _n = other._n;
            _taps = _alloc.allocate(_n);
            std::uninitialized_copy_n(other._taps, _n, _taps);
        }
        return *this;
    }
    constexpr fir_filter& operator=(fir_filter&& other){
        if(this != other){
            _n = other._n;
            _taps = std::move(other._taps);
            other._taps = nullptr;
        }
        return *this;
    }


    constexpr size_type getNumTaps() const noexcept{return _n;}
    constexpr pointer getTaps() const noexcept{return _taps;}

    constexpr void design(const size_type& sps, const value_type& beta, const size_type& span) 
    requires(__F == filter_type::Root_Raised_Cosine){
        if(sps < 1.0) throw std::logic_error("Samples per symbol cannot be less than 1");
        if(span <= 0.0) throw std::logic_error("Span must be greater than 0 symbols");
        if(beta <= 0.0 || beta > 1.0) throw std::logic_error("Beta must be in range (0, 1]");

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
        for(size_type i = 0; i < n; ++i){
            _taps[i] /= power;
        }
    }

    constexpr void design(const size_type& order, const value_type& fc)
    requires(__F == filter_type::Lowpass_Ideal){
        _taps = _alloc.allocate(order);
        _n = order;
        difference_type n_max = order % 2 ? (order - 1)/2 : order/2;
        if(order % 2) _taps[n_max] = fc / PI;
        for(difference_type n = 0; n < n_max; ++n){
            difference_type m = n - n_max;
            _taps[n] = std::sin(m*fc)/m/PI;
            _taps[order - 1 - n] = _taps[n];
        }
    }

    constexpr void design(const size_type& order, const value_type& fc)
    requires(__F == filter_type::Highpass_Ideal){
        _taps = _alloc.allocate(order);
        _n = order;
        difference_type n_max = order % 2 ? (order - 1)/2 : order/2;
        if(order % 2) _taps[n_max] = 1.0 - fc / PI;
        for(difference_type n = 0; n < n_max; ++n){
            difference_type m = n - n_max;
            _taps[n] = -std::sin(m*fc)/m/PI;
            _taps[order - 1 - n] = _taps[n];
        }
    }

    constexpr void design(const size_type& order, const value_type& fc1, const value_type& fc2)
    requires(__F == filter_type::Bandpass_Ideal){
        _taps = _alloc.allocate(order);
        _n = order;
        difference_type n_max = order % 2 ? (order - 1)/2 : order/2;
        if(order % 2) _taps[n_max] = (fc2 - fc1)/PI;
        for(difference_type n = 0; n < n_max; ++n){
            difference_type m = n - n_max;
            _taps[n] = (std::sin(m*fc2) - std::sin(m*fc1))/m/PI;
            _taps[order - 1 - n] = _taps[n];
        }
    }

    constexpr void design(const size_type& order, const value_type& fc1, const value_type& fc2)
    requires(__F == filter_type::Bandstop_Ideal){
        _taps = _alloc.allocate(order);
        _n = order;
        difference_type n_max = order % 2 ? (order - 1)/2 : order/2;
        if(order % 2) _taps[n_max] = 1.0 + (fc1 - fc2)/PI;
        for(difference_type n = 0; n < n_max; ++n){
            difference_type m = n - n_max;
            _taps[n] = (std::sin(m*fc1) - std::sin(m*fc2))/m/PI;
            _taps[order - 1 - n] = _taps[n];
        }
    }

    constexpr void design(const size_type& order, const pointer samples)
    requires(__F == filter_type::Frequency_Sampled){
        _n = order;
        _taps = _alloc.allocate(_n);
        size_type M = (_n - 1.0)/2.0;
        for(size_type n = 0; n < order; ++n){
            value_type tap = samples[0];
            value_type x = 2.0*PI*(static_cast<value_type>(n) - M)/_n;
            for(size_type k = 1; k <= (_n - 1)/2; ++k) tap += 2.0*std::cos(x*k)*samples[k];
            _taps[n] = tap/_n;
        }
    }

    constexpr void synthesize()
    requires(__W == window_type::None){} // also a no-op, but is here for standards, so calling design() synthesize() always works

    constexpr void synthesize()
    requires(__W == window_type::Rectangular){} // does nothing lol, but needs to be here, probably optimized to no-op anyways

    constexpr void synthesize()
    requires(__W == window_type::Triangular){
        value_type midpoint = static_cast<value_type>(_n - 1)/2.0;
        for(size_type n = 0; n < _n; ++n){
            _taps[n] *= 1.0 - std::abs((n - midpoint)/midpoint);
        }
    }

    constexpr void synthesize()
    requires(__W == window_type::Von_Hann){
        for(size_type n = 0; n < _n; ++n){
            _taps[n] *= 0.5*(1.0 - std::cos((2.0*PI*n)/(_n-1)));
        }
    }

    constexpr void synthesize()
    requires(__W == window_type::Hamming){
        for(size_type n = 0; n < _n; ++n){
            _taps[n] *= 0.54 - 0.46*std::cos((2.0*PI*n)/(_n-1));
        }
    }

    constexpr void synthesize()
    requires(__W == window_type::Blackman_Harris){
        for(size_type n = 0; n < _n; ++n){
            _taps[n] *= 0.35875 - 0.48829*std::cos(2.0*PI*n/_n) + 0.14128*std::cos(4.0*PI*n/_n) - 0.01168*std::cos(6.0*PI*n/_n);
        }
    }

    constexpr void synthesize(const value_type& sigma = 0.4)
    requires(__W == window_type::Gaussian){
        for(size_type n = 0; n < _n; ++n){
            _taps[n] *= std::exp(-0.5*std::pow((n - _n/2.0)/sigma/_n/2.0, 2.0));
        }
    }

    constexpr void synthesize(const value_type& alpha = 0.4)
    requires(__W == window_type::Tukey){
        size_type m = _n - 1;
        for(size_type n = 0; n < _n; ++n){
            if(n < alpha*m/2.0){
                _taps[n] *= 0.5*(1 + std::cos(PI*(2.0*n)/alpha/m - 1.0));
            } else if(n > (_n - 1.0)*(1.0 - alpha/2.0)){
                _taps[n] *= 0.5*(1 + std::cos(PI*(2.0*n)/alpha/m - 1.0));
            }
        }
    }

    constexpr std::pair<pointer, pointer> response(const size_type& points){
        std::pair<pointer, pointer> out;
        out.first = _alloc.allocate(points);
        out.second = _alloc.allocate(points);
        for(size_type idx = 0; idx < points; ++idx){
            __T lambda = idx*PI/points;
            std::complex<__T> work(0.0, 0.0);
            for(size_type tap_idx = 0; tap_idx < _n; ++tap_idx){
                work += _taps[tap_idx] * std::complex<__T>(std::cos(tap_idx*lambda), -1.0*std::sin(tap_idx*lambda));
            }
            out.first[idx] = lambda;
            out.second[idx] = 20.0 * std::log10(std::norm(work));
        }
        return out;
    }

    constexpr pointer filter(pointer input, const size_type& M){
        pointer out = _alloc.allocate(M);
        for(size_type m = 0; m < M; ++m){
            value_type work = value_type();
            for(size_type n = 0; n < _n; ++n){
                if(m >= n) work += _taps[n] * input[m - n];
            }
            out[m] = work;
        }
        return out;
    }

    constexpr iterator begin() noexcept{return iterator(&_taps[0]);}
    constexpr const_iterator begin() const noexcept{return const_iterator(&_taps[0]);}
    constexpr const_iterator cbegin() const noexcept{return const_iterator(&_taps[0]);}
    
    constexpr iterator end() noexcept{return iterator(&_taps[_n]);}
    constexpr const_iterator end() const noexcept{return const_iterator(&_taps[_n]);}
    constexpr const_iterator cend() const noexcept{return const_iterator(&_taps[_n]);}
    
    constexpr reverse_iterator rbegin() noexcept{return iterator(&_taps[_n - 1]);}
    constexpr const_reverse_iterator rbegin() const noexcept{return const_reverse_iterator(&_taps[_n - 1]);}
    constexpr const_reverse_iterator crbegin() const noexcept{return const_reverse_iterator(&_taps[_n - 1]);}
    
    constexpr reverse_iterator rend() noexcept{return iterator(&_taps[-1]);}
    constexpr const_reverse_iterator rend() const noexcept{return const_reverse_iterator(&_taps[-1]);}
    constexpr const_reverse_iterator crend() const noexcept{return const_reverse_iterator(&_taps[-1]);}
    
};
#endif