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

#define DEBUG

#include <complex>
#include <memory>
#include <algorithm>
#include <ranges>
#include <vector>
#include <numbers>
#include <numeric>
#include <cmath>
#include <bitset>

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
        out_t out =  operate();
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
                accum += _taps[j] * std::polar(1.0, -lambda * j);
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
        double* temp = new double[_prot_size];
        std::uninitialized_copy(rg.begin(), rg.end(), temp);
        update_taps(temp, _prot_size);
        _state = INTERP;
    }

    /**
     * @brief Update the resampler with a range compatible container
     * @param rate resampling rate
     * @param num_filters total arms in the bank
     * @param rg range containing the prototype
     */
    template<class R>
    requires (std::ranges::input_range<R> && std::convertible_to<std::ranges::range_reference_t<R>, double>)
    void update(const double& rate, const size_t& num_filters, R&& rg){
        std::destroy_n(_bank, _n_filts);
        std::destroy_n(_deriv_bank, _n_filts);
        delete[] _bank;
        delete[] _deriv_bank;
        _rate = rate;
        _n_filts = num_filters;
        _tau = 0.0;
        _theta = 0.0;
        _mu = 0.0;
        _curr = std::ceil(_n_filts / 2);
        x = 0.0;
        dx = 0.0;
        _bank = new bank_t[_n_filts];
        _deriv_bank = new bank_t[_n_filts];
        _prot_size = std::ranges::distance(rg);
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
        _prot_size = n;
        _taps_per_filter = std::ceil(_prot_size / _n_filts);
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
        double differentiator[3] = {-1.0, 0.0, 1.0};
        double* diff_filt = new double[_prot_size];
        double power = 0.0;
        for(size_t i = 0; i < _prot_size - 2; ++i){
            double accum = 0;
            for(size_t j = 0; j < 3; ++j){
                accum += differentiator[j] * taps[i + j];
            }
            diff_filt[i] = accum;
            power += std::abs(accum);
        }
        diff_filt[_prot_size - 1] = 0.0;
        auto f = [this, power](double& tap){tap *= _n_filts / power;};
        std::for_each(diff_filt, diff_filt + _prot_size, f);

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
        reset();
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
        _state = INTERP;
        for(size_t i = 0; i < _n_filts; ++i){
            _bank[i].reset();
            _deriv_bank[i].reset();
        }
    }

    eftc<double, out_t>* get_bank() const{return _bank;}
    eftc<double, out_t>* get_deriv_bank() const{return _deriv_bank;}
    double get_rate() const{return _rate;}
    size_t get_num_filters() const{return _n_filts;}
    size_t get_taps_per_arm() const{return _taps_per_filter;}
    size_t get_sample_delay() const{return _taps_per_filter;}

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

    /**
     * @brief Set the arm to use on the next operation
     * @param arm 
     */
    void set_next_arm(size_t arm){_curr = arm;}
};

/**
 * @brief Rectangular modulator is a generic modulator for rectangular constellations 
 */
class rectangular_modulator{
private:
    constellation* _constel; // rectangular constellation to use
    size_t _bits_per_symbol; // number of bits per symbol, determines input buffer length
    size_t _sps; // samples per symbol
    size_t _n_filts;
    size_t _n_points; // number of points in the constellation
    double _rolloff; // matched filter rolloff alpha

    std::vector<double> _prototype; // prototype filter used
    resampler<std::complex<double>> _pfb; // polyphase resampler to use
public:
    rectangular_modulator() = delete;

    /**
     * @brief Construct a new modulator
     * @param constellation the constellation to use
     * @param buffer_size the total number of complex symbols to hold in an internal buffer
     * @param sps samples per symbol to use
     */
    rectangular_modulator(constellation* constellation, size_t sps, size_t num_filters, double rolloff)
        : _constel(constellation), _bits_per_symbol(_constel->get_bps()), _sps(sps), _n_filts(num_filters)
        ,_n_points(constellation->get_size()), _rolloff(rolloff)
        , _prototype(root_nyquist(num_filters, num_filters, 1.0, _rolloff, 8 * _sps * num_filters))
        , _pfb(_sps, num_filters, _prototype){
    }

    constellation* get_constellation() const{return _constel;}
    size_t get_bps() const{return _bits_per_symbol;}
    size_t get_sps() const{return _sps;}
    std::vector<double> get_prototype() const{return _prototype;}
    resampler<std::complex<double>> get_resampler() const{return _pfb;}

    void set_constellation(constellation* constel){_constel = constel; update_internals();}
    void set_sps(size_t sps){_sps = sps; update_internals();}
    void set_num_filters(size_t n){_n_filts = n; update_internals();}
    void set_rolloff(double alpha){_rolloff = alpha; update_internals();}

    /**
     * @brief Generate sps symbols per bps bits
     * @param input an integer representing an index into the constellation
     * @param output buffer to hold complex symbols, must be at least length sps
     * @return int number of symbols placed into the output buffer
     */
    int operate(size_t input, std::complex<double>* output){
        if(input >= _n_points) throw std::out_of_range("MERM operate: input sample is out of range");
        auto point = _constel->get_point(input);
        int n = _pfb.operate(point, output);
        return n;
    }

    void reset(){
        _pfb.reset();
    }

    void update_internals(){
        _prototype.clear();
        _prototype = root_nyquist(_n_filts, _n_filts, 1.0, _rolloff, 8 * _sps * _n_filts);
        _pfb.update(_sps, _n_filts, _prototype);
        reset();
    }

};

#ifdef DEBUG
/**
 * @brief Contains arrays of all internal loop variables
 */
struct DEBUG_INTERFACE{
    size_t n = 500;
    size_t curr = 0;

    /* AGC VARIABLES */
    double AGC_RATE;
    double AGC_REF;
    double AGC_GAIN;
    double AGC_MAX_GAIN;

    double* AGC_GAIN_HIST = new double[n];

    /* FLL VARIABLES */
    double FLL_DAMPING, FLL_ALPHA, FLL_BETA, FLL_MIN_FREQ, FLL_MAX_FREQ;
    size_t FLL_PROTO_SIZE;

    double* FLL_ERR_HIST = new double[n];
    double* FLL_PHASE_HIST = new double[n];
    double* FLL_FREQ_HIST = new double[n];

    /* PLL VARIABLES */
    double PLL_DAMPING, PLL_ALPHA, PLL_BETA, PLL_MIN_FREQ, PLL_MAX_FREQ;
    
    double* PLL_ERR_HIST = new double[n];
    double* PLL_PHASE_HIST = new double[n];
    double* PLL_FREQ_HIST = new double[n];
};
#endif

/**
 * @brief Firefighter is a full digital demodulation class incorporating frequency, phase, and symbol recovery loops
 * 
 */
class firefighter{
private:
    #ifdef DEBUG
    DEBUG_INTERFACE _debug;
    #endif

    size_t _sps; // samples per symbol of the input complex stream
    size_t _n_filts; // number of filters to use for matched filtering
    constellation* _constel; // the constellation to use
    double _loop_bw; // internal loop bandwidth for control
    double _rolloff; // matched filter rolloff factor alpha

    static constexpr double PI = std::numbers::pi;

    /* AUTOMATIC GAIN CONTROL VARIABLES */
    double _agc_rate;
    double _agc_ref;
    double _agc_gain;
    double _agc_max_gain;
    /************************************/

    /* FREQUENCY RECOVERY VARIABLES */
    double _fll_phase, _fll_freq;
    double _fll_max_freq, _fll_min_freq;
    double _fll_damping;
    double _fll_alpha, _fll_beta;
    size_t _fll_size;

    double _fll_err;

    eftc<std::complex<double>, std::complex<double>> _fll_lowerband_filter;
    eftc<std::complex<double>, std::complex<double>> _fll_upperband_filter;
    /********************************/

    /* PREAMBLE ESTIMATION VARIABLES */
    enum PEV_STATE{
        COLLECT, // Collect samples for correlation
        ONE_MORE, // PEV has found a point of correlation and needs one more sample
        OPERATE,
        FINISHED,
        WAIT
    } _pev_state;

    std::complex<double>* _pev_sync_words;
    size_t _pev_sync_size;
    size_t* _pev_preamble;
    size_t _pev_preamble_size;
    size_t _pev_mark_delay, _pev_old_mark_delay;
    double _pev_threshold, _pev_old_threshold;
    double _pev_pfa;
    size_t _pev_corr_idx;

    size_t _pev_corr_est;
    double _pev_time_est;
    double _pev_phase_est;

    rectangular_modulator* _pev_mod;

    std::complex<double>* _pev_corr_hist;
    double* _pev_corr_mag_hist;
    eftc<std::complex<double>, std::complex<double>> _pev_filter;
    /********************************/

    /* PHASE RECOVERY VARIABLES */
    double _pll_phase, _pll_freq;
    double _pll_max_freq, _pll_min_freq;
    double _pll_damping;
    double _pll_alpha, _pll_beta;

    double _pll_err;
    std::complex<double>* _pll_output_buffer;

    std::vector<double> _pll_prototype;
    resampler<std::complex<double>> _pll_pfb;
    /***************************/

    /* DEBUG INTERFACE */
    double* _pll_err_debug = new double[500];
    double* _pll_phase_debug = new double[500];
    double* _pll_freq_debug = new double[500];
    /*******************/

    double sinc(double x){
        if(x > -1e-6 && x < 1e-6) return 1.0;
        else return std::sin(PI * x) / (PI * x);
    }
    
    size_t* uint8_to_points(const std::bitset<8>& val){
        const size_t n = _constel->get_bps();
        size_t k = 8 / n;
        size_t* output = new size_t[k];
        std::bitset<8> mask = 0xFF;
        mask >>= (8 - n);
        for(size_t i = 0; i < k; ++i){
            std::bitset<8> x = (val & (mask << (i * n))) >> (i * n);
            output[i] = x.to_ullong();
            // std::cout<<"i: "<<i<<", val: "<<val<<", mask: "<<(mask << (i * n))<<", x: "<<x<<", output: "<<output[i]<<"\n";
        }
        return output;
    }

public:
    firefighter() = delete;

    /**
     * @brief Construct a new recovery object
     * @param constel the constellation used to modulate the data
     * @param sps the input stream samples per second
     * @param n the number of filters to use for matched filter decimation
     * @param loop_bandwidth the bandwidth of the internal control loop
     */
    firefighter(constellation* constel, const size_t& sps, const size_t& n, double loop_bandwidth, double rolloff)
        : _sps(sps), _n_filts(n), _constel(constel), _loop_bw(loop_bandwidth), _rolloff(rolloff)
        , _pll_prototype(root_nyquist(_n_filts, _n_filts * _sps, 1.0, _rolloff, 8 * _n_filts * _sps))
        , _pll_pfb(1.0 / _sps, _n_filts, _pll_prototype){
        update_internals();
    }

    size_t get_sps() const{return _sps;}
    size_t get_num_filters() const{return _n_filts;}
    double get_bandwidth() const{return _loop_bw;}
    constellation* get_constellation() const{return _constel;}

    size_t get_fll_size() const{return _fll_size;}
    eftc<std::complex<double>, std::complex<double>> get_fll_lower_band() const{return _fll_lowerband_filter;}
    eftc<std::complex<double>, std::complex<double>> get_fll_upper_band() const{return _fll_upperband_filter;}
    size_t get_sample_delay() const{return _pll_pfb.get_sample_delay();}

    size_t* get_preamble_points() const{return _pev_preamble;}
    size_t get_preamble_size() const{return _pev_preamble_size;}
    std::complex<double>* get_sync_word() const{return _pev_sync_words;}
    size_t get_sync_size() const{return _pev_sync_size;}
    double get_pev_threshold() const{return _pev_threshold;}

    void set_sps(size_t sps){_sps = sps; update_internals();}
    void set_num_filters(size_t n){_n_filts = n; update_internals();}
    void set_bandwidth(double bandwidth){_loop_bw = bandwidth; update_internals();}
    void set_constellation(constellation* constel){_constel = constel; update_internals();}

    void set_fll_size(size_t size){_fll_size = size; update_fll();}

    std::vector<double> get_pll_prototype() const{return _pll_prototype;}

    DEBUG_INTERFACE* debug(){return &_debug;}

    /**
     * @brief Run the AGC on an input complex sample
     * @param input 
     * @return std::complex<double> 
     */
    std::complex<double> agc(std::complex<double> input){
        std::complex<double> agc_output = input * _agc_gain;
        _agc_gain += _agc_rate * (_agc_ref - std::sqrt(std::pow(agc_output.real(), 2.0) + std::pow(agc_output.imag(), 2.0)));
        if(_agc_gain > _agc_max_gain){
            _agc_gain = _agc_max_gain;
        }
        #ifdef DEBUG
        _debug.AGC_GAIN_HIST[_debug.curr] = _agc_gain;
        #endif
        return agc_output;
    }

    /**
     * @brief Run the FLL on an input complex sample
     * @param input 
     * @return std::complex<double> 
     */
    std::complex<double> fll(std::complex<double> input){
        std::complex<double> fll_nco = std::polar(1.0,  _fll_phase);
        std::complex<double> fll_output = input * fll_nco;

        std::complex<double> out_lower = _fll_lowerband_filter.filter(fll_output);
        std::complex<double> out_upper = _fll_upperband_filter.filter(fll_output);

        _fll_err = std::norm(out_upper) - std::norm(out_lower);
        _fll_freq += _fll_beta * _fll_err;
        _fll_phase += _fll_freq;

        if(_fll_phase > 2.0 * PI) _fll_phase = std::fmod(_fll_phase, 2.0 * PI);
        if(_fll_phase < -2.0 * PI) _fll_phase = std::fmod(_fll_phase, -2.0 * PI);
        _fll_freq = _fll_freq > _fll_max_freq ? _fll_max_freq : _fll_freq;
        _fll_freq = _fll_freq < _fll_min_freq ? _fll_min_freq : _fll_freq;

        #ifdef DEBUG
        _debug.FLL_ERR_HIST[_debug.curr] = _fll_err;
        _debug.FLL_FREQ_HIST[_debug.curr] = _fll_freq;
        _debug.FLL_PHASE_HIST[_debug.curr] = _fll_phase;
        #endif
        return fll_output;
    }

    /**
     * @brief Run the PLL on an input complex sample, which acts as a pure decimator
     * @param input 
     * @return std::complex<double> 
     */
    std::complex<double> pll(std::complex<double> input, int& n){
        n = _pll_pfb.operate(input, _pll_output_buffer);
        std::complex<double> pll_output = 0.0;
        if(n != 0){
            if(_pev_state == FINISHED){
                // pev has a phase estimate for us, so use it
                pll_output = _pll_output_buffer[0] * std::polar(1.0, -_pev_phase_est);
                _pev_state = WAIT;
            } 
            else{
                std::complex<double> pll_nco = std::polar(1.0, -_pll_phase);
            pll_output = *_pll_output_buffer * pll_nco;
            }

            _pll_err = _constel->phase_error_detector(pll_output);

            _pll_freq += _pll_beta * _pll_err;
            _pll_phase += _pll_freq + _pll_alpha * _pll_err;

            if(_pll_phase > 2.0 * PI) _pll_phase = std::fmod(_pll_phase, 2.0 * PI);
            if(_pll_phase < -2.0 * PI) _pll_phase = std::fmod(_pll_phase, -2.0 * PI);
            _pll_freq = _pll_freq > _pll_max_freq ? _pll_max_freq : _pll_freq;
            _pll_freq = _pll_freq < _pll_min_freq ? _pll_min_freq : _pll_freq;
        }
        #ifdef DEBUG
        _debug.PLL_ERR_HIST[_debug.curr] = _pll_err;
        _debug.PLL_PHASE_HIST[_debug.curr] = _pll_phase;
        _debug.PLL_FREQ_HIST[_debug.curr] = _pll_freq;
        #endif
        return pll_output;
    }

    /**
     * @brief Runs the PEV on an input, which may modify the internal state and pass timing and phase estimates to a section
     * @param input 
     * @return true if the PEV has found an estimate
     */
    bool pev(std::complex<double> input){
        if(_pev_state == FINISHED || _pev_state == WAIT) return true;

        _pev_filter.feed_sample(input);
        std::complex<double> pev_corr = _pev_filter.operate();
        _pev_filter.increment();

        double pev_corr_mag = std::norm(pev_corr);

        _pev_corr_hist[_pev_corr_idx] = pev_corr;
        _pev_corr_mag_hist[_pev_corr_idx] = pev_corr_mag;

        if(pev_corr_mag <= _pev_threshold && _pev_state == COLLECT){
            // We havent reached a good correlation so we should continue sampling
            _pev_corr_idx++;
            _pev_state = COLLECT;
            return false;
        }

        if(_pev_state == COLLECT){
            // correlation has reached its maximum, but we need 1 more sample
            _pev_state = ONE_MORE;
            return false;
        }
        if(_pev_state == ONE_MORE){
            // _pev_corr_idx++;
            _pev_state = OPERATE;
        }
        // We reach here only if correlation > threshold and state = ONE_MORE || OPERATE
        _pev_time_est = (_pev_corr_mag_hist[_pev_corr_idx - 1] + 2.0 * _pev_corr_mag_hist[_pev_corr_idx] + 3.0 * _pev_corr_mag_hist[_pev_corr_idx + 1])
                        / (_pev_corr_mag_hist[_pev_corr_idx - 1] + _pev_corr_mag_hist[_pev_corr_idx] + _pev_corr_mag_hist[_pev_corr_idx + 1]) - 2.0;
        _pev_phase_est = std::atan2(_pev_corr_hist[_pev_corr_idx].imag(), _pev_corr_hist[_pev_corr_idx].real());
        _pev_corr_est = _pev_corr_idx + _pev_mark_delay;

        _pev_state = FINISHED;

        _pev_corr_idx += static_cast<size_t>(_sps + 0.5);
        return true;
    }

    /**
     * @brief Operate the symbol recovery on one sample.
     * @param sample complex sample from input stream
     * @param output buffer having enough space for the output
     */
    int operate(const std::complex<double>& sample, std::complex<double>* output){
        int n = 0; // number of samples written to the output in this cycle (decimation notifier)
        std::complex<double> agc_out = agc(sample);
        std::complex<double> fll_out = fll(agc_out);

        // preamble detected should be used only if the state is correct
        bool detected = pev(fll_out);

        // if(detected){
        //     // the pll will take the estimate given
        //     std::complex<double> pll_out = pll(fll_out, n);
        //     output[0] = pll_out;
        // }
        std::complex<double> pll_out = pll(fll_out, n);
            output[0] = pll_out;

        #ifdef DEBUG
        _debug.curr = (_debug.curr + 1) % _debug.n;
        #endif
        return n;
    }

    void reset(){
        _fll_upperband_filter.reset();
        _fll_lowerband_filter.reset();
        _pll_pfb.reset();
        _pev_filter.reset();
        _pev_state = COLLECT;
        _pev_corr_idx = 0;
    }

    void set_pll_alpha(double alpha){_pll_alpha = alpha;}
    void set_pll_beta(double beta){_pll_beta = beta;}
    void set_max_freq(double max){_pll_max_freq = max;}
    void set_min_freq(double min){_pll_min_freq = min;}

private:
    void update_agc(){
        _agc_rate = 1e-4;
        _agc_ref = 1.0;
        _agc_gain = 1.0;
        _agc_max_gain = 1e3;
    }

    void update_fll(){
        _fll_damping = std::sqrt(2.0) / 2.0;
        _fll_alpha = 0.0;
        _fll_beta = 4.0 * _loop_bw / static_cast<double>(_sps);
        _fll_max_freq = 2.0 * PI / static_cast<double>(_sps);
        _fll_min_freq = -2.0 * PI / static_cast<double>(_sps);
        _fll_phase = 0.0;
        _fll_freq = 0.0;
        _fll_size = 55;
        _fll_err = 0.0;
        set_fll_filter();
    }

    void update_preamble(){
        // Update the reference to the modulator used, and setup the sync words
        _pev_mod = new rectangular_modulator(_constel, _sps, _n_filts, _rolloff);
        const size_t preamble_len = 8;
        const size_t k = 8 / _constel->get_bps();
        const std::bitset<8> preamble[preamble_len] = {0xAA, 0x00, 0xFE, 0x2E, 0x22, 0xFC, 0x05, 0x70};
        _pev_preamble = new size_t[k * preamble_len];
        _pev_preamble_size = k * preamble_len;
        for(size_t i = 0; i < preamble_len; ++i){
            size_t* points = uint8_to_points(preamble[i]);
            for(size_t j = 0; j < k; ++j){
                _pev_preamble[i * k + j] = points[j];
            }
        }

        // Setup the modulator sync word
        _pev_sync_size = _sps * _pev_preamble_size;
        _pev_sync_words = new std::complex<double>[_pev_sync_size];
        for(size_t i = 0; i < _pev_preamble_size; ++i){
            _pev_mod->operate(_pev_preamble[i], _pev_sync_words + i * _sps);
        }
        for(size_t i = 0; i < _pev_sync_size; ++i){
            _pev_sync_words[i] = std::conj(_pev_sync_words[i]);
        }
        std::reverse(_pev_sync_words, _pev_sync_words + _pev_sync_size);
        _pev_filter.update_taps(_pev_sync_words, _pev_sync_size);
        _pev_corr_hist = new std::complex<double>[_pev_sync_size];
        _pev_corr_mag_hist = new double[_pev_sync_size];

        _pev_mark_delay = 1;
        _pev_old_mark_delay = 1;
        _pev_threshold = 0.9;
        _pev_old_threshold = 0.9;
        _pev_pfa = -std::log(1.0 - _pev_threshold);
        _pev_corr_idx = 0;

        _pev_state = COLLECT;

        double accum = 0.0;
        for(size_t i = 0; i < _pev_sync_size; ++i){
            accum += std::abs(_pev_sync_words[i] * std::conj(_pev_sync_words[i]));
        }
        _pev_threshold = _pev_threshold * std::pow(accum, 2.0);
    }

    void update_pll(){
        _pll_prototype = root_nyquist(_n_filts, _n_filts * _sps, 1.0, _rolloff, 8 * _n_filts * _sps);
        _pll_pfb.update(1.0 / _sps, _n_filts, _pll_prototype);

        _pll_phase = 0.0;
        _pll_freq = 0.0;
        _pll_max_freq = 0.25;
        _pll_min_freq = -1.0 * _pll_max_freq;
        _pll_damping = std::sqrt(2.0) / 2.0;
        _pll_alpha = 4.0 * _pll_damping * _loop_bw / (1.0 + 2.0 * _pll_damping * _loop_bw + std::pow(_loop_bw, 2.0));
        _pll_beta = 4.0 * std::pow(_loop_bw, 2.0) / (1.0 + 2.0 * _pll_damping * _loop_bw + std::pow(_loop_bw, 2.0));
        _pll_err = 0.0;
        _pll_output_buffer = new std::complex<double>[_sps];
        
        set_pll_filter();
    }

    void set_fll_filter(){
        long long int M = std::rint(static_cast<double>(_fll_size) / static_cast<double>(_sps));
        double power = 0.0;

        double* temp = new double[_fll_size];
        for(size_t i = 0; i < _fll_size; ++i){
            const double k = static_cast<double>(-M) + (static_cast<double>(i) * 2.0 / static_cast<double>(_sps));
            const double pos = _rolloff * k;
            double tap = sinc(pos + 0.5) + sinc(pos - 0.5);
            power += std::pow(tap, 2.0);
            temp[i] = tap;
        }

        std::complex<double>* upper = new std::complex<double>[_fll_size];
        std::complex<double>* lower = new std::complex<double>[_fll_size];
        long long int N = (_fll_size - 1) / 2;
        for(size_t i = 0; i < _fll_size; ++i){
            double tap = temp[i] / power;
            double k = (static_cast<double>(i) - static_cast<double>(N)) * 0.5 / static_cast<double>(_sps);
            size_t idx = _fll_size - i - 1;
            lower[idx] = std::polar(tap, -2.0 * PI * (1.0 + _rolloff) * k);
            upper[idx] = std::conj(lower[idx]);
        }
        _fll_lowerband_filter.update_taps(lower, _fll_size);
        _fll_upperband_filter.update_taps(upper, _fll_size);
    }

    void set_pll_filter(){
    }

    #ifdef DEBUG
    void setup_debug(){
        _debug.AGC_GAIN = _agc_gain;
        _debug.AGC_MAX_GAIN = _agc_max_gain;
        _debug.AGC_RATE = _agc_rate;
        _debug.AGC_REF = _agc_ref;

        _debug.FLL_DAMPING = _fll_damping;
        _debug.FLL_ALPHA = _fll_alpha;
        _debug.FLL_BETA = _fll_beta;
        _debug.FLL_PROTO_SIZE = _fll_size;
        _debug.FLL_MIN_FREQ = _fll_min_freq;
        _debug.FLL_MAX_FREQ = _fll_max_freq;

        _debug.PLL_DAMPING = _pll_damping;
        _debug.PLL_ALPHA = _pll_alpha;
        _debug.PLL_BETA = _pll_beta;
        _debug.PLL_MIN_FREQ = _pll_min_freq;
        _debug.PLL_MAX_FREQ = _pll_max_freq;

        std::uninitialized_default_construct_n(_debug.AGC_GAIN_HIST, _debug.n);
        std::uninitialized_default_construct_n(_debug.FLL_ERR_HIST, _debug.n);
        std::uninitialized_default_construct_n(_debug.FLL_FREQ_HIST, _debug.n);
        std::uninitialized_default_construct_n(_debug.FLL_PHASE_HIST, _debug.n);
        std::uninitialized_default_construct_n(_debug.PLL_ERR_HIST, _debug.n);
        std::uninitialized_default_construct_n(_debug.PLL_FREQ_HIST, _debug.n);
        std::uninitialized_default_construct_n(_debug.PLL_PHASE_HIST, _debug.n);
    }
    #endif

    /**
     * @brief Update all internal parameters of the recovery because of updated parameters
     */
    void update_internals(){
        update_fll(); // initialize fll internal variables and update fll filters
        update_preamble();
        update_agc(); // initialize the AGC section
        update_pll(); // initialize the pll and its polyphase bank
        reset(); // reset internal states for filters
        #ifdef DEBUG
        setup_debug();
        #endif
    }
};

#endif