/**
 * @file modem.hpp
 * @author Caleb Brown (cwbrown@mines.edu)
 * @brief Main Header file for all digital modulation recovery, Real time operating.
 * @version 0.1
 * @date 2025-10-13
 * @copyright Copyright (c) 2025
 */
#ifndef MODEM_HPP
#define MODEM_HPP

#include <complex>
#include <memory>
#include <algorithm>
#include <ranges>
#include <vector>
#include <complex>
#include <numbers>
#include <numeric>
#include <cmath>

/**
 * @brief Frequency response struct, containing the sample frequencies and response magnitude in dB
 */
struct FREQ_RESP{
    double* frequencies;
    double* magnitudes;
    size_t n;
};

/**
 * @brief Memory efficient FIR filter class
 * @tparam tap_t the type used for taps
 * @tparam out_t the type of samples used
 */
template<typename tap_t, typename out_t>
class eftc{
private:
    size_t _N; // number of taps
    tap_t* _taps; // taps array
    size_t _curr; // current index into the history
    out_t* _history; // filter history

    static constexpr double PI = std::numbers::pi;
public:
    /**
     * @brief Construct a null filter
     */
    eftc(): _N(0), _taps(nullptr), _curr(0), _history(nullptr){}

    /**
     * @brief Construct a new FIR filter with taps from a pointer
     * @param taps pointer to the new taps
     * @param N number of taps in the prototype
     */
    eftc(tap_t* taps, size_t N): _N(N), _taps(new tap_t[_N]), _curr(0), _history(new out_t[_N]){
        std::uninitialized_copy_n(taps, _N, _taps);
        std::uninitialized_default_construct_n(_history, _N);
    }

    /**
     * @brief Construct an FIR filter with an iterator 
     * @param start an iterator to the start of a container
     * @param end  an iterator to the end of a container
     */
    template<class InputIt>
    eftc(InputIt start, InputIt end)
        : _N(std::distance(start, end)), _taps(new tap_t[_N]), _curr(0), _history(new out_t[_N]){
        std::uninitialized_copy(start, end, _taps);
        std::uninitialized_default_construct_n(_history, _N);
    }

    /**
     * @brief Initialize the filter with a range compatible container
     * @param range the range containing the filter taps
     */
    template<class R>
    requires (std::ranges::input_range<R> && std::convertible_to<std::ranges::range_reference_t<R>, tap_t>)
    eftc(R&& range)
        : _N(std::ranges::distance(range)), _taps(new tap_t[_N]), _curr(0), _history(new out_t[_N]){
        std::uninitialized_copy(range.begin(), range.end(), _taps);
        std::uninitialized_default_construct_n(_history, _N);
    }

    ~eftc(){}

    size_t get_num_taps() const{return _N;}
    tap_t* get_taps() const{return _taps;}

    /**
     * @brief Filter a single sample once
     * @param sample 
     * @return std::complex<double> 
     */
    out_t filter(const out_t& sample){
        _history[_curr] = sample;
        out_t accum(0.0);
        size_t idx = _curr;
        for(size_t i = 0; i < _N; ++i){
            accum += _taps[i] * _history[idx];
            idx = idx == 0 ? _N - 1 : idx - 1;
        }
        _curr = (_curr + 1ULL) % _N;
        return accum;
    }

    out_t stepback_filter(){
        out_t accum(0.0);
        size_t idx = _curr;
        for(size_t i = 0; i < _N; ++i){
            idx = idx == 0 ? _N - 1 : idx - 1;
            accum += _taps[i] * _history[idx];
        }
        return accum;
    }

    /**
     * @brief Reset the internal state of the filter
     */
    void reset(){
        std::fill_n(_history, _N, out_t{});
        _curr = 0;
    }

    /**
     * @brief Update the filter with new taps
     * @param taps pointer to the new taps
     * @param n number of taps
     */
    void update_taps(tap_t* taps, size_t n){
        std::destroy_n(_taps, _N);
        std::destroy_n(_history, _N);
        delete[] _taps;
        delete[] _history;
        _N = n;
        _curr = 0;
        _taps = new tap_t[_N];
        _history = new out_t[_N];
        std::uninitialized_copy_n(taps, _N, _taps);
        std::uninitialized_default_construct_n(_history, _N);
    }

    /**
     * @brief Update the filter with new taps from a range
     * @param range the range containing the new taps
     */
    template<class R>
    requires (std::ranges::input_range<R> && std::convertible_to<std::ranges::range_reference_t<R>, tap_t>)
    void update_taps(R&& range){
        std::destroy_n(_taps, _N);
        std::destroy_n(_history, _N);
        delete[] _taps;
        delete[] _history;
        _N = std::ranges::distance(range);
        _curr = 0;
        _taps = new tap_t[_N];
        _history = new out_t[_N];
        std::uninitialized_copy(range.begin(), range.end(), _taps);
        std::uninitialized_default_construct_n(_history, _N);
    }

    /**
     * @brief Calculate the frequency domain response of the filter
     * 
     * @param points 
     * @return FREQ_RESP 
     */
    FREQ_RESP response(const size_t& points){
        FREQ_RESP output;
        output.frequencies = new double[points];
        output.magnitudes = new double[points];
        output.n = points;
        for(size_t i = 0; i < points; ++i){
            double lambda = i * PI / points;
            std::complex<double> accum(0.0, 0.0);
            for(size_t j = 0; j < _N; ++j){
                accum += std::polar(_taps[j], -lambda * j);
            }
            output.frequencies[i] = lambda;
            output.magnitudes[i] = 10.0 * std::log10(std::norm(accum));
        }
        return output;
    }
};

/**
 * @brief Polyphase arbitrary resampling filter
 * @tparam out_t the type of samples to filter
 */
template<typename out_t>
class resampler{
private:
    size_t _interp; // interpolation rate
    size_t _decim; // decimation rate
    size_t _n; // number of filters in the bank
    size_t _prot_size; // number of taps in the prototype

    double _accum; // fractional accumulation of the resample
    double _rate; // input rate of the resample
    double _frac_rate; // fractional rate of the resample
    size_t _taps_per_filter; // taps per polyphase filter arm bank
    size_t _curr; // current bank arm

    eftc<double, out_t>* _bank;
    eftc<double, out_t>* _deriv_bank;

    // out_t* _history;
    // size_t _hist_curr;

    static constexpr double PI = std::numbers::pi;
public:
    /**
     * @brief Construct a null polyphase resampling filter
     */
    resampler(): _bank(nullptr), _deriv_bank(nullptr){}

    /**
     * @brief Construct a resampler with taps from a pointer
     * @param rate resampling rate
     * @param taps pointer to an array of taps
     * @param n size of the prototype
     * @param num_filters number of filters in the polyphase bank
     */
    resampler(const double& rate, double* taps, const size_t& n, const size_t& num_filters)
        : _interp(num_filters), _decim(std::floor(_interp / rate)), _n(num_filters)
        , _prot_size(n), _accum(0.0), _rate(rate)
        , _frac_rate(_interp / rate - _decim), _taps_per_filter(std::ceil(_prot_size / _interp))
        , _curr(std::ceil(_n / 2)) 
        , _bank(new eftc<double, out_t>[_n])
        , _deriv_bank(new eftc<double, out_t>[_n]){
        // Create a temporary array containing the prototype padded with zeros
        double* temp = new double[_interp * _taps_per_filter];
        std::uninitialized_default_construct_n(temp, _interp * _taps_per_filter);
        std::uninitialized_copy_n(taps, n, temp);

        // Create our polyphase bank
        for(size_t i = 0; i < _n; ++i){
            // A temporary array to hold this phases taps
            double* temp_sub = new double[_taps_per_filter];
            for(size_t j = 0; j < _taps_per_filter; ++j){
                temp_sub[j] = temp[i + j * _interp];
            }
            _bank[i].update_taps(temp_sub, _taps_per_filter);
        }

        // Create the derivative of our prototype
        double differentiator[2] = {-1.0, 1.0};
        double* _diff_filt = new double[_prot_size];
        for(size_t i = 0; i < _prot_size - 1; ++i){
            double accum = 0;
            for(size_t j = 0; j < 2; ++j){
                accum += differentiator[j] * taps[i + j]; // weird notation but it does work on valid iterators
            }
            _diff_filt[i] = accum;
        }
        _diff_filt[_prot_size - 1] = 0.0;

        // Now create the derivative bank
        double* deriv_temp = new double[_interp * _taps_per_filter];
        std::uninitialized_default_construct_n(temp, _interp * _taps_per_filter);
        std::uninitialized_copy_n(_diff_filt, _prot_size, deriv_temp);
        for(size_t i = 0; i < _n; ++i){
            // A temporary array to hold this phases taps
            double* temp_sub = new double[_taps_per_filter];
            for(size_t j = 0; j < _taps_per_filter; ++j){
                temp_sub[j] = deriv_temp[i + j * _interp];
            }
            _deriv_bank[i].update_taps(temp_sub, _taps_per_filter);
        }
    }

    /**
     * @brief Construct a resampler with taps from the range [start, end]
     * @param rate resampling rate
     * @param start iterator to start of range
     * @param end iterator to end of range
     * @param num_filters number of filters in the polyphase bank
     */
    template<class InputIt>
    resampler(const double& rate, InputIt start, InputIt end, const size_t& num_filters)
        : _interp(num_filters), _decim(std::floor(_interp / rate)), _n(num_filters)
        , _prot_size(std::distance(start, end)), _accum(0.0), _rate(rate)
        , _frac_rate(_interp / rate - _decim), _taps_per_filter(std::ceil(_prot_size / _interp))
        , _curr(std::ceil(_n / 2)) 
        , _bank(new eftc<double, out_t>[_n])
        , _deriv_bank(new eftc<double, out_t>[_n]){
        // Create a temporary array containing the prototype padded with zeros
        double* temp = new double[_interp * _taps_per_filter];
        std::uninitialized_default_construct_n(temp, _interp * _taps_per_filter);
        std::uninitialized_copy(start, end, temp);

        // Create our polyphase bank
        for(size_t i = 0; i < _n; ++i){
            // A temporary array to hold this phases taps
            double* temp_sub = new double[_taps_per_filter];
            for(size_t j = 0; j < _taps_per_filter; ++j){
                temp_sub[j] = temp[i + j * _interp];
            }
            _bank[i].update_taps(temp_sub, _taps_per_filter);
        }

        // Create the derivative of our prototype
        double differentiator[2] = {-1.0, 1.0};
        double* _diff_filt = new double[_prot_size];
        for(size_t i = 0; i < _prot_size - 1; ++i){
            double accum = 0;
            for(size_t j = 0; j < 2; ++j){
                accum += differentiator[j] * *(start + i + j); // weird notation but it does work on valid iterators
            }
            _diff_filt[i] = accum;
        }
        _diff_filt[_prot_size - 1] = 0.0;

        // Now create the derivative bank
        double* deriv_temp = new double[_interp * _taps_per_filter];
        std::uninitialized_default_construct_n(temp, _interp * _taps_per_filter);
        std::uninitialized_copy_n(_diff_filt, _prot_size, deriv_temp);
        for(size_t i = 0; i < _n; ++i){
            // A temporary array to hold this phases taps
            double* temp_sub = new double[_taps_per_filter];
            for(size_t j = 0; j < _taps_per_filter; ++j){
                temp_sub[j] = deriv_temp[i + j * _interp];
            }
            _deriv_bank[i].update_taps(temp_sub, _taps_per_filter);
        }
    }

    /**
     * @brief Construct a resampler with taps from a compatible range
     * @param rate resampling rate
     * @param range range containing prototype taps
     * @param num_filters number of filters in the polyphase bank
     */
    template<class R>
    requires (std::ranges::input_range<R> && std::convertible_to<std::ranges::range_reference_t<R>, double>)
    resampler(const double& rate, R&& range, const size_t& num_filters)
        : _interp(num_filters), _decim(std::floor(_interp / rate)), _n(num_filters)
        , _prot_size(std::ranges::distance(range)), _accum(0.0), _rate(rate)
        , _frac_rate(_interp / rate - _decim), _taps_per_filter(std::ceil(_prot_size / _interp))
        , _curr(std::ceil(_n / 2)) 
        , _bank(new eftc<double, out_t>[_n])
        , _deriv_bank(new eftc<double, out_t>[_n]){
        // Create a temporary array containing the prototype padded with zeros
        double* temp = new double[_interp * _taps_per_filter];
        std::uninitialized_default_construct_n(temp, _interp * _taps_per_filter);
        std::uninitialized_copy(range.begin(), range.end(), temp);

        // Create our polyphase bank
        for(size_t i = 0; i < _n; ++i){
            // A temporary array to hold this phases taps
            double* temp_sub = new double[_taps_per_filter];
            for(size_t j = 0; j < _taps_per_filter; ++j){
                temp_sub[j] = temp[i + j * _interp];
            }
            _bank[i].update_taps(temp_sub, _taps_per_filter);
        }

        // Create the derivative of our prototype
        double differentiator[2] = {-1.0, 1.0};
        double* _diff_filt = new double[_prot_size];
        for(size_t i = 0; i < _prot_size - 1; ++i){
            double accum = 0;
            for(size_t j = 0; j < 2; ++j){
                accum += differentiator[j] * range[i + j]; // weird notation but it does work on valid iterators
            }
            _diff_filt[i] = accum;
        }
        _diff_filt[_prot_size - 1] = 0.0;

        // Now create the derivative bank
        double* deriv_temp = new double[_interp * _taps_per_filter];
        std::uninitialized_default_construct_n(temp, _interp * _taps_per_filter);
        std::uninitialized_copy_n(_diff_filt, _prot_size, deriv_temp);
        for(size_t i = 0; i < _n; ++i){
            // A temporary array to hold this phases taps
            double* temp_sub = new double[_taps_per_filter];
            for(size_t j = 0; j < _taps_per_filter; ++j){
                temp_sub[j] = deriv_temp[i + j * _interp];
            }
            _deriv_bank[i].update_taps(temp_sub, _taps_per_filter);
        }
    }

    eftc<double, out_t>* get_bank() const{return _bank;}
    eftc<double, out_t>* get_deriv_bank() const{return _deriv_bank;}
    size_t get_interpolation() const{return _interp;}
    size_t get_decimation() const{return _decim;}
    size_t get_num_filters() const{return _n;}
    size_t get_taps_per_arm() const{return _taps_per_filter;}

    /**
     * @brief Operate the resampler on a sample value
     * @param sample 
     * @param output a buffer having enough space for 
     * @return bool true when a sample is put into the output buffer
     */
    bool operate(const out_t& sample, out_t* output){
        for(size_t i = 0; i < _n; ++i){
            _bank[i].filter(sample);
            _deriv_bank[i].filter(sample);
        }

        // while(_curr < _n){
            
        // }
    }

};

#endif