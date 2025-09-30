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

/**
 * @brief Returns the taps for a root nyquist filter.
 * @param gain total gain
 * @param fs sampling frequency
 * @param sr symbol rate (symbols/s)
 * @param alpha excess bandwidth/rolloff 
 * @param n number of taps
 * @return std::vector<double> 
 */
std::vector<double> root_nyquist(const double& gain, const double& fs, const double& sr, const double& alpha, const size_t& n){
    if(gain <= 0.0) throw std::domain_error("make_nyquist: gain cannot be less than 0");
    if(alpha <= 0.0 || alpha > 1.0) throw std::domain_error("make_nyquist: rolloff must be in range (0, 1]");

    const size_t n_taps = n | 1;
    std::vector<double> taps(n_taps, 0.0);
    double power = 0.0;
    const double sps = fs/sr;
    const double PI = std::numbers::pi;

    for(size_t i = 0; i < n_taps; ++i){
        double num, denom;
        double xidx = (double)i - (double)n_taps/2.0;
        double x1 = PI*xidx/(double)sps;
        double x2 = 4.0*alpha*xidx/(double)sps;
        double x3 = std::pow(x2, 2.0) - 1.0;

        if(std::abs(x3) >= 1.0e-5){
            if(i != n_taps/2){
                num = std::cos((1.0 + alpha)*x1) 
                    + std::sin((1.0 - alpha)*x1)/(4.0*alpha*xidx/(double)sps);
            }
            else{
                num = std::cos((1.0 + alpha)*x1)
                    + (1.0 - alpha)*PI/(4.0*alpha);
            }
            denom = x3*PI;
        }
        else{
            if(alpha == 1.0){
                taps[i] = -1.0;
                power += taps[i];
                continue;
            }
            x3 = (1.0 - alpha)*x1;
            x2 = (1.0 + alpha)*x1;
            num = std::sin(x2)*(1.0 + alpha)*PI
                - std::cos(x3)*((1.0 - alpha)*PI*(double)sps)/(4.0*alpha*xidx)
                + std::sin(x3)*std::pow((double)sps, 2.0)/(4.0*alpha*std::pow(xidx, 2.0));
            denom = -32.0*PI*std::pow(alpha, 2.0)*xidx/(double)sps;
        }
        taps[i] = 4.0*alpha*num/denom;
        power += taps[i];
    }
    for(size_t i = 0; i < n_taps; ++i){
        taps[i] = taps[i]*gain/power;
    }
    return taps;
}

/**
 * @brief Container class for taps of an FIR filter.
 * 
 */
class fir_taps{
private:
    std::vector<double> _taps;
    size_t _n;
    std::vector<std::complex<double>> _history;
    size_t _curr = 0;
public:
    fir_taps(const std::vector<double>& taps, const size_t& n): _taps(taps), _n(n), _history(n, 0.0){}
    std::complex<double> filter1(const std::complex<double>& sample){
        std::complex<double> accum(0.0, 0.0);
        _history[_curr] = sample;
        size_t idx = _curr;
        for(size_t i = 0; i < _n; ++i){
            accum += _taps[i]*_history[idx];
            if(idx == 0) idx = _history.size() - 1;
            else idx -= 1;
        }
        _curr = (_curr + 1)%_history.size();
        return accum;
    }
    void update_taps(const std::vector<double>& taps){_taps = taps; _n = taps.size();}

    void reset(){
        _history.assign(_n, 0.0);
        _curr = 0;
    }

    std::vector<double> taps() const{return _taps;}
    size_t size() const{return _n;}
};

/**
 * @brief polyphase_filter_bank is an internal class of symbol recovery that peforms as a polyphase filter bank. It is made with a prototype, 
 * usually a root nyquist filter.
 */
class polyphase_filter_bank{
private:
    size_t _n; // number of arms of the bank
    size_t _taps_per_filter;
    std::vector<fir_taps> _filters;
    std::vector<fir_taps> _deriv_filters;
    std::vector<double> _prototype;
    std::vector<double> _deriv_prototype;

    void set_taps(){
        std::vector<double> prot_filled(_prototype);
        while(prot_filled.size() < _n * _taps_per_filter) prot_filled.push_back(0.0);
        for(size_t i = 0; i < _n; ++i){
            std::vector<double> tmp(_taps_per_filter, 0.0);
            for(size_t j = 0; j < _taps_per_filter; ++j){
                tmp[j] = prot_filled[i + j*_n];
            }
            _filters[i].update_taps(tmp);
        }

        std::vector<double> deriv_filled(_deriv_prototype);
        while(deriv_filled.size() < _n * _taps_per_filter) deriv_filled.push_back(0.0);
        for(size_t i = 0; i < _n; ++i){
            std::vector<double> tmp(_taps_per_filter, 0.0);
            for(size_t j = 0; j < _taps_per_filter; ++j){
                tmp[j] = deriv_filled[i + j*_n];
            }
            _deriv_filters[i].update_taps(tmp);
        }
    };

    void make_derivative(){
        _deriv_prototype.reserve(_prototype.size());
        _deriv_prototype.push_back(0.0);
        std::vector<double> diff{-1.0, 0.0, 1.0};
        double power = 0.0;
        for(size_t i = 0; i < _prototype.size() - 2; ++i){
            double accum = 0.0;
            for(size_t j = 0; j < 3; ++j){
                accum += diff[j]*_prototype[i + j];
            }
            _deriv_prototype.push_back(accum);
            power += std::abs(accum);
        }
        _deriv_prototype.push_back(0.0);
        for(size_t i = 0; i < _deriv_prototype.size(); ++i){
            _deriv_prototype[i] *= _n/power;
        }
    };
public:
    polyphase_filter_bank(): _n(0), _taps_per_filter(0){}
    constexpr polyphase_filter_bank(const std::vector<double>& prototype, const size_t& arms){
        update(prototype, arms);
    }

    size_t get_size() const{return _n;}
    size_t get_taps_per_arm() const{return _taps_per_filter;}
    fir_taps get_arm(const size_t& idx) const{return _filters.at(idx);}
    fir_taps get_deriv_arm(const size_t& idx) const{return _deriv_filters.at(idx);}
    std::vector<double> get_prototype() const{return _prototype;}
    std::vector<double> get_deriv_prototype() const{return _deriv_prototype;}

    std::complex<double> filter(const std::complex<double>& sample, const size_t& arm){
        return _filters[arm].filter1(sample);
    }
    std::complex<double> deriv_filter(const std::complex<double>& sample, const size_t& arm){
        return _deriv_filters[arm].filter1(sample);
    }

    constexpr void update(const std::vector<double>& prototype, const size_t& arms){
        _n = arms;
        _prototype = prototype;
        make_derivative();

        // Set up our tap vectors
        _taps_per_filter = std::ceil(prototype.size()/_n);
        _filters.reserve(_n);
        _deriv_filters.reserve(_n);
        fir_taps tmp_taps(std::vector<double>(1, 0.0), 1);
        for(size_t i = 0; i < _n; ++i){
            _filters.emplace_back(tmp_taps);
            _deriv_filters.emplace_back(tmp_taps);
        }

        set_taps();
    }

    void clear(){
        _filters.clear();
        _deriv_filters.clear();
        _prototype.clear();
        _deriv_prototype.clear();
    }
};

class symbol_recovery{
private:
    // *** Storage Parameters ***
    size_t _sps;
    double _bandwidth;
    double _filter_bandwidth;
    size_t _n_filters;
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
    size_t _curr; // The current arm of the bank
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
        std::vector<double> prototype = root_nyquist(_n_filters, _n_filters, 1.0/_sps, _filter_bandwidth, 8*_sps*_n_filters);
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
    void operate(std::vector<std::complex<double>> samples){
        size_t i = 0;
        size_t j = 0;
        while(i < samples.size()){
            _curr = std::floor(_phase);

            // Wrap around our filter bank to the proper arm

            auto v = _bank.filter(samples[j], _curr);
            _phase = _phase + _deviation;

            auto dv = _bank.deriv_filter(samples[j], _curr);
            double real_err = v.real() * dv.real();
            double imag_err = v.imag() * dv.imag();
            _error = (real_err + imag_err) / 2.0;

            for(size_t s = 0; s < _sps; ++s){
                _deviation += _ki * _error;
                _phase += _deviation + _kp * _error;
            }

            _deviation = 0.5 * (std::abs(_deviation + _max_deviation) - std::abs(_deviation - _max_deviation));

            i++;
            j += _sps;
        }
    }
};


#endif