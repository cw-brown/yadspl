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
#include <numbers>
#include <numeric>
#include <cmath>

#include "constellations.hpp"

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
     * @brief Push one sample into the internal history
     * @param sample 
     */
    void feed_sample(const out_t& sample){
        _history[_curr] = sample;
    }

    /**
     * @brief Run the filter on the current history of samples
     * @return out_t 
     */
    out_t operate(){
        out_t accum(0.0);
        size_t idx = _curr;
        for(size_t i = 0; i < _N; ++i){
            accum+= _taps[i] * _history[idx];
            idx = idx == 0 ? _N - 1 : idx - 1;
        }
        return accum;
    }

    /**
     * @brief Increment the internal history counter
     */
    void increment(){
        _curr = (_curr + 1) % _N;
    }

    /**
     * @brief Filter a single sample once
     * @param sample 
     * @return std::complex<double> 
     */
    out_t filter(const out_t& sample){
        feed_sample(sample);
        auto out =  operate();
        increment();
        return out;
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
public:
    using bank_t = eftc<double, out_t>;
private:
    enum SAMP_STATE{
        BOUNDARY,
        INTERP
    } _state;

    double _rate; // resampling rate
    size_t _n_filts; // number of filters in the bank

    double _tau; // phase accumulation
    double _theta; // soft filter index, theta=tau*n
    double _mu; // fractional interpolation rate
    size_t _curr; // current arm of the bank

    size_t _taps_per_filter; // number of taps per bank arm filter
    size_t _prot_size; // number of taps in the prototype

    eftc<double, out_t>* _bank;
    eftc<double, out_t>* _deriv_bank;

    out_t x;
    out_t dx;

    static constexpr double PI = std::numbers::pi;
public:
    /**
     * @brief Construct a null polyphase resampling filter
     */
    resampler(): _bank(nullptr), _deriv_bank(nullptr){}

    /**
     * @brief Construct a resampler with a pointer to a prototype
     * @param rate resampling rate
     * @param num_filters total arms in the bank
     * @param taps pointer to prototype filter
     * @param n number of taps in the prototype
     */
    resampler(const double& rate, const size_t& num_filters, double* taps, size_t n)
        : _rate(rate), _n_filts(num_filters), _tau(0.0), _theta(0.0), _mu(0.0), _curr(std::ceil(_n_filts / 2))
        , _bank(new bank_t[_n_filts]), _deriv_bank(new bank_t[_n_filts]), x(0.0), dx(0.0){
        _prot_size = n;
        _taps_per_filter = std::ceil(_prot_size / num_filters);
        update_taps(taps, n);
        _state = INTERP;
    }

    /**
     * @brief Construct a resampler with a range given by iterators
     * @param rate resampling rate
     * @param num_filters total arms in the bank
     * @param start iterator pointing to the start of the prototype
     * @param end iterator pointing to the end of the prototype
     */
    template<class InputIt>
    resampler(const double& rate, const size_t& num_filters, InputIt start, InputIt end)
        : _rate(rate), _n_filts(num_filters), _tau(0.0), _theta(0.0), _mu(0.0), _curr(std::ceil(_n_filts / 2))
        , _bank(new bank_t[_n_filts]), _deriv_bank(new bank_t[_n_filts]), x(0.0), dx(0.0){
        _prot_size = std::distance(start, end);
        _taps_per_filter = std::ceil(_prot_size / num_filters);
        double* temp = new double[_prot_size];
        std::uninitialized_copy(start, end, temp);
        update_taps(temp, _prot_size);
        _state = INTERP;
    }

    /**
     * @brief Construct a resampler with a range compatible container
     * @param rate resampling rate
     * @param num_filters total arms in the bank
     * @param rg range containing the prototype
     */
    template<class R>
    requires (std::ranges::input_range<R> && std::convertible_to<std::ranges::range_reference_t<R>, double>)
    resampler(const double& rate, const size_t& num_filters, R&& rg)
        : _rate(rate), _n_filts(num_filters), _tau(0.0), _theta(0.0), _mu(0.0), _curr(std::ceil(_n_filts / 2))
        , _bank(new bank_t[_n_filts]), _deriv_bank(new bank_t[_n_filts]), x(0.0), dx(0.0){
        _prot_size = std::ranges::distance(rg);
        _taps_per_filter = std::ceil(_prot_size / num_filters);
        double* temp = new double[_prot_size];
        std::uninitialized_copy(rg.begin(), rg.end(), temp);
        update_taps(temp, _prot_size);
        _state = INTERP;
    }

    /**
     * @brief Update the internal polyphase bank of the filter with a prototype
     * @param taps pointer to an FIR prototype taps
     * @param n size of the prototype
     */
    void update_taps(double* taps, const size_t& n){
        // Allocate a temporary buffer to store a zero padded array of the prototype
        double* temp = new double[_n_filts * _taps_per_filter];
        std::uninitialized_default_construct_n(temp, _n_filts * _taps_per_filter);
        std::uninitialized_copy_n(taps, n, temp);

        // Create our polyphase bank
        for(size_t i = 0; i < _n_filts; ++i){
            // A temporary array to hold this phases taps
            double* temp_sub = new double[_taps_per_filter];
            for(size_t j = 0; j < _taps_per_filter; ++j){
                temp_sub[j] = temp[i + j * _n_filts];
            }
            _bank[i].update_taps(temp_sub, _taps_per_filter);
        }

        // Create a derivative of our prototype
        double differentiator[2] = {-1.0, 1.0};
        double* diff_filt = new double[_prot_size];
        for(size_t i = 0; i < _prot_size - 1; ++i){
            double accum = 0;
            for(size_t j = 0; j < 2; ++j){
                accum += differentiator[j] * taps[i + j];
            }
            diff_filt[i] = accum;
        }
        diff_filt[_prot_size - 1] = 0.0;

        // Create the derivative polyphase bank
        double* deriv_temp = new double[_n_filts * _taps_per_filter];
        std::uninitialized_default_construct_n(temp, _n_filts * _taps_per_filter);
        std::uninitialized_copy_n(diff_filt, _prot_size, deriv_temp);
        for(size_t i = 0; i < _n_filts; ++i){
            // A temporary array to hold this phases taps
            double* temp_sub = new double[_taps_per_filter];
            for(size_t j = 0; j < _taps_per_filter; ++j){
                temp_sub[j] = deriv_temp[i + j * _n_filts];
            }
            _deriv_bank[i].update_taps(temp_sub, _taps_per_filter);
        }
    }

    /**
     * @brief Change the resampling rate
     * 
     * @param rate 
     */
    void set_rate(const double& rate){
        _rate = rate;
    }

    /**
     * @brief Reset the internal state of the resampler
     */
    void reset(){
        _tau = 0.0;
        _theta = 0.0;
        _mu = 0.0;
        _curr = std::ceil(_n_filts / 2);
        x = 0.0;
        dx = 0.0;
    }

    eftc<double, out_t>* get_bank() const{return _bank;}
    eftc<double, out_t>* get_deriv_bank() const{return _deriv_bank;}
    double get_rate() const{return _rate;}
    size_t get_num_filters() const{return _n_filts;}
    size_t get_taps_per_arm() const{return _taps_per_filter;}

    /**
     * @brief Operate the resampler on a sample value
     * @param sample 
     * @param output a buffer having enough space for all outputs
     * @return int number of samples written to the output buffer
     */
    int operate(const out_t& sample, out_t* output){
        // Add the sample internally to all arms
        for(size_t i = 0; i < _n_filts; ++i){
            _bank[i].feed_sample(sample);
            _deriv_bank[i].feed_sample(sample);
        }

        size_t n = 0;

        while(_curr < _n_filts){
            switch(_state){
            case BOUNDARY:
                dx = _deriv_bank[_curr].operate();
                output[n++] = (1.0 - _mu) * x + _mu * dx;
                _tau += 1.0 / _rate;
                _theta = _tau * static_cast<double>(_n_filts);
                _curr = std::floor(_theta);
                _mu = _theta - _curr;
                _state = INTERP;
                break;
            case INTERP:
                x = _bank[_curr].operate();
                if(_curr == _n_filts - 1){
                    _state = BOUNDARY;
                    _curr = _n_filts;
                }
                else{
                    x = _bank[_curr + 1].operate();
                    output[n++] = (1.0 - _mu) * x + _mu * dx;
                    _tau += 1.0 / _rate;
                    _theta = _tau * static_cast<double>(_n_filts);
                    _curr = std::floor(_theta);
                    _mu = _theta - _curr;
                }
                break;
            default:
                return -1;
            }
        }

        _tau -= 1.0f;
        _theta -= static_cast<double>(_n_filts);
        _curr -= _n_filts;

        for(size_t i = 0; i < _n_filts; ++i){
            _bank[i].increment();
            _deriv_bank[i].increment();
        }
        return n;
    }

};

/**
 * @brief Memory Efficient Rectangular Modulator (MERM) is a generic modulator for rectangular constellations 
 */
class merm{
private:
    constellation* _constel; // rectangular constellation to use
    size_t _bits_per_symbol; // number of bits per symbol, determines input buffer length
    size_t _sps; // samples per symbol

    std::vector<double> _prototype; // prototype filter used
    resampler<std::complex<double>> _pfb; // polyphase resampler to use
public:
    merm() = delete;

    /**
     * @brief Construct a new modulator
     * @param constellation the constellation to use
     * @param buffer_size the total number of complex symbols to hold in an internal buffer
     * @param sps samples per symbol to use
     */
    merm(constellation* constellation, const size_t& buffer_size, size_t sps, size_t num_filters)
        : _constel(constellation), _bits_per_symbol(_constel->get_bps()), _sps(sps)
        , _prototype(root_nyquist(num_filters, num_filters, 1.0, 0.35, 8 * _sps * num_filters))
        , _pfb(_sps, num_filters, _prototype){
    }

    constellation* get_constellation() const{return _constel;}
    size_t get_bps() const{return _bits_per_symbol;}
    size_t get_sps() const{return _sps;}
    std::vector<double> get_prototype() const{return _prototype;}
    resampler<std::complex<double>> get_resampler() const{return _pfb;}

    /**
     * @brief Generate sps symbols per bps bits
     * @param input an integer representing an index into the constellation
     * @param output buffer to hold complex symbols, must be at least length sps
     * @return int number of symbols placed into the output buffer
     */
    int operate(unsigned int input, std::complex<double>* output){
        if(input >= _bits_per_symbol) throw std::out_of_range("MERM operate: input sample is out of range");
        auto point = _constel->get_point(input);
        int n = _pfb.operate(point, output);
        return n;
    }

};

#endif