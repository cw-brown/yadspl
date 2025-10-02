/**
 * @file symbol_rec.hpp
 * @author Caleb Brown
 * @brief Full symbol recovery for complex symbols.
 * @version 0.1
 * @date 2025-09-18
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#ifndef SYMBOL_REC_H
#define SYMBOL_REC_H

#include <vector>
#include <complex>
#include <numbers>
#include <numeric>
#include <cmath>

#include "polyphase.hpp"

/**
 * @brief 
 * 
 */
class carrier_recovery{
private:
    size_t _sps; // samples per symbol
    double _alpha; // input filter roll-off
    double _bandwidth; // loop bandwidth
    size_t _n; // prototype filter size

    double _phase;
    double _freq;
    double _max_freq;
    double _min_freq;
    double _damping;
    double _alpha;
    double _beta;

    fir_taps _lower_band;
    fir_taps _upper_band;

    double sinc(double x){
        return x == 0.0 ? 1.0 : std::sin(PI * x) / x;
    }

    void update_filter(){
        int M = std::round(_n / _sps);
        double pow = 0.0;

        std::vector<double> baseband;
        baseband.reserve(_n);
        const double half = 2.0 / _sps;
        for(size_t i = 0; i < _n; ++i){
            const double k = -M + i * half;
            const double pos = _alpha * k;
            const double tap = sinc(pos - 0.5) + sinc(pos + 0.5);
            pow += std::pow(tap, 2.0);
            baseband.push_back(tap);
        }

        std::vector<std::complex<double>> _upper(_n);
        std::vector<std::complex<double>> _lower(_n);

        long int N = (baseband.size() - 1) / 2;
        for(size_t i = 0; i < _n; ++i){
            const double tap = baseband[i] / pow;
            const double k = (static_cast<int>(i) - N) * 0.5 / _sps;
            size_t idx = _n - i - 1;
            _lower[idx] = tap * std::exp(-2.0 * PI * (1.0 + _alpha) * k);
            _upper[idx] = std::conj(_lower[idx - i - 1]);
        }

    }

    static constexpr double PI = std::numbers::pi;
public:
    carrier_recovery(const size_t& sps, const double& roll_off, const double& loop_bw, const size_t& filter_size)
        : _sps(sps), _alpha(roll_off), _bandwidth(loop_bw), _n(filter_size){
        // Set up control variables and gains
        _phase = 0;
        _freq = 0;
        _min_freq = -4.0 * PI / _sps;
        _max_freq = 4.0 * PI / _sps;
        _damping = std::sqrt(2.0) / 2.0;
        _alpha = 0.0;
        _beta = 8.0 * PI * _bandwidth / _sps;
        update_filter();
    }

    size_t sps() const{return _sps;}
    double damp() const{return _damping;}
    double max_freq() const{return _max_freq;}
    double min_freq() const{return _min_freq;}
    double alpha() const{return _alpha;}
    double beta() const{return _beta;}
    double phase() const{return _phase;}
    double frequency() const{return _freq;}

};


class symbol_recovery{
private:
    // *** Storage Parameters ***
    size_t _sps;
    double _bandwidth;
    double _filter_bandwidth;
    int _n_filters;
    polyphase_filter_bank _bank;
    constellation _constel;

    // *** Loop Parameters ***
    double _phase; // The current phase offset (also the filter number)
    double _max_deviation; // The maximum rate at which the loop can move from 0
    double _deviation; // The current deviation from 0
    double _damp; // Dampening coefficient
    double _kp; // Proportional gain
    double _ki; // Integral gain
    double _error; // Total error signal
    int _curr; // The current arm of the bank
public:
    /**
     * @brief Symbol recovery will receive an input sample stream and output decoded symbols (bytes) according to a provided constellation.
     * The object creates a filter bank that is matched to the transmitter's filter. Symbol Recovery is a very memory intensive object
     * and there should only ever be one at a time. It owns all its components.
     * @param sps The samples-per-symbol of the input stream
     * @param loop_bandwidth Bandwidth and gain of the control loop
     * @param num_filters The total number of filters within the recovery block
     * @param max_deviation Maximum rate at which the loop can deviate from 0
     * @param constellation A constellation object that encodes the data
     * @param filter_bandwidth The bandwidth of the matched filter
     */
    symbol_recovery(const size_t& sps, const double& loop_bandwidth, const size_t& num_filters, 
                    const double& max_deviation, const constellation& constellation, const double& filter_bandwidth)
        : _sps(sps), _bandwidth(loop_bandwidth), _filter_bandwidth(filter_bandwidth), 
        _n_filters(num_filters), _constel(constellation), _phase(num_filters/2.0),
        _max_deviation(max_deviation), _deviation(0.0), _damp(0.0), _kp(0.0), _ki(0.0), _error(0.0), _curr(0)
    {
        // Generate the polyphase prototype
        std::vector<double> prototype = root_nyquist(_n_filters, _n_filters*_sps, 1.0, _filter_bandwidth, 8*_sps*_n_filters);
        _bank = polyphase_filter_bank(prototype, _n_filters);

        // Update internal parameters of the control loop
        update_internals();
        _curr = std::floor(_phase);
    }

    size_t get_sps() const{return _sps;}
    double get_loop_bandwidth() const{return _bandwidth;}
    size_t get_num_filters() const{return _n_filters;}
    double get_phase() const{return _phase;}
    double get_max_deviation() const{return _max_deviation;}
    double get_deviation() const{return _deviation;}
    double get_dampening() const{return _damp;}
    double get_p_gain() const{return _kp;}
    double get_i_gain() const{return _ki;}
    double get_error() const{return _error;}
    size_t get_current_arm() const{return _curr;}
    polyphase_filter_bank get_bank() const{return _bank;}

    void set_sps(const size_t& sps){
        _sps = sps;
        update_filter();
    }
    void set_filter_bandwidth(const double& bandwidth){
        _filter_bandwidth = bandwidth;
        update_filter();
    }
    void set_loop_bandwidth(const double& bandwidth){
        _bandwidth = bandwidth;
        update_internals();
    }
    void set_num_filters(const size_t& n){
        _n_filters = n;
        update_filter();
        update_internals();
    }
    void set_max_deviation(const double& deviation){
        _max_deviation = deviation;
    }
    void set_constellation(const constellation& constel){
        _constel = constel;
    }

    void update_filter(){
        std::vector<double> prototype = root_nyquist(_n_filters, _n_filters, 1.0/_sps, _filter_bandwidth, 8*_sps*_n_filters);
        _bank.clear();
        _bank.update(prototype, prototype.size());
    }
    void update_internals(){
        _damp = 2.0 * _n_filters;
        const double denom = 1.0 + 2.0 * _damp * _bandwidth + std::pow(_bandwidth, 2.0);
        _kp = (4.0 * _damp * _bandwidth) / denom;
        _ki = (4.0 * std::pow(_bandwidth, 2.0)) / denom;
    }

    /**
     * @brief Runs the internal loop on a collection of samples
     * 
     * @param samples 
     */
    std::vector<std::complex<double>> operate(std::vector<std::complex<double>> samples){
        std::vector<std::complex<double>> output;
        for(size_t i = 0; i < samples.size(); i+= _sps){
            _curr = std::floor(_phase);

            if(_curr >= _n_filters){
                _phase -= _n_filters;
                _curr %= _n_filters;
            }
            while(_curr < 0){
                _phase += _n_filters;
                _curr += _n_filters;
            }



        }
        return output;
    }

};


#endif