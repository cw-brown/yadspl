#ifndef MODEM2_HPP
#define MODEM2_HPP

#include <cmath>
#include <numeric>
#include <numbers>
#include <bitset>

#include "constellations.hpp"


template<typename tap_t, typename out_t>
class fir{
private:
    size_t _n;
    tap_t* _taps;
    size_t _curr;
    out_t* _history;
public:
    constexpr fir(): _n(0), _taps(nullptr), _curr(0), _history(nullptr){}
    constexpr fir(tap_t* taps, const size_t& n): _n(n), _taps(new tap_t[_n]), _curr(0), _history(new out_t[_n]){
        std::uninitialized_copy_n(taps, _n, _taps);
        std::uninitialized_default_construct_n(_history, _n);
    }
    ~fir(){delete[] _taps; delete[] _history;}

    constexpr size_t num_taps() const noexcept{return _n;}
    constexpr tap_t* taps() const noexcept{return _taps;}

    constexpr void set_taps(tap_t* taps, const size_t& n) noexcept{
        if(_taps != nullptr) delete[] _taps;
        if(_history != nullptr) delete[] _history;
        _n = n; _curr = 0;
        _taps = new tap_t[n];
        _history = new out_t[n];
        std::uninitialized_copy_n(taps, _n, _taps);
        std::uninitialized_default_construct_n(_history, _n);
    }
    
    constexpr void insert(const out_t& sample){
        _history[_curr] = sample;
    }
    constexpr out_t operate(){
        out_t accum(0.0);
        size_t idx = _curr;
        for(size_t i = 0; i < _n; ++i){
            accum += _taps[i] * _history[idx];
            idx = idx == 0 ? _n - 1 : idx - 1;
        }
        return accum;
    }
    constexpr void increment(){
        _curr = (_curr + 1) % _n;
    }

    constexpr void filter(out_t input, out_t* output){
        insert(input);
        output[0] = operate();
        increment();
    }
    constexpr void filter_n(out_t* input, const size_t& n_in, out_t* output, size_t& n_written){
        for(size_t i = 0; i < n_in; ++i){
            insert(input[i]);
            output[i] = operate();
            increment();
        }
        n_written = n_in;
    }
    constexpr out_t filter_norm(out_t input){
        insert(input);
        out_t output = operate();
        increment();
        return output;
    }

    constexpr void reset(){
        std::fill_n(_history, _n, out_t{});
        _curr = 0;
    }
};

template<typename tap_t, typename out_t>
class fir_pfb{
private:
    size_t _prot_size; // total number of filter coefficients
    size_t _arm_size; // total number of coefficients in each arm
    size_t _num_filters;

    out_t* _history;
    size_t _curr;
    tap_t** _taps;
public:
    constexpr fir_pfb():_history(nullptr), _taps(nullptr){}
    constexpr fir_pfb(size_t num_filters, tap_t* taps, size_t n){
        _curr = 0;
        _taps = nullptr;
        _history = nullptr;
        update(num_filters, taps, n);
    }
    ~fir_pfb(){
        delete[] _history;
        for(size_t i = 0; i < _num_filters; ++i) delete[] _taps[i];
        delete[] _taps;
    }

    constexpr size_t prototype_size() const noexcept{return _prot_size;}
    constexpr size_t arm_size() const noexcept{return _arm_size;}
    constexpr size_t num_filters() const noexcept{return _num_filters;}
    constexpr tap_t* arm(size_t idx) const noexcept{return _taps[idx];}

    constexpr void update(size_t num_filters, tap_t* taps, size_t n){
        _num_filters = num_filters;
        _prot_size = n;
        _arm_size = std::ceil(_prot_size / _num_filters);

        if(_taps != nullptr){
            for(size_t i = 0; i < _num_filters; ++i) delete[] _taps[i];
            delete[] _taps;
        }
        _taps = new tap_t*[_num_filters];
        for(size_t i = 0; i < _num_filters; ++i) _taps[i] = new tap_t[_arm_size];
        if(_history != nullptr) delete[] _history;
        _history = new out_t[_arm_size];
        std::uninitialized_default_construct_n(_history, _arm_size);

        tap_t* temp = new tap_t[_num_filters * _arm_size];
        std::uninitialized_default_construct_n(temp, _num_filters * _arm_size);
        std::uninitialized_copy_n(taps, _prot_size, temp);

        for(size_t i = 0; i < _num_filters; ++i){
            for(size_t j = 0; j < _arm_size; ++j){
                _taps[i][j] = temp[i + j * _num_filters];
            }
        }
        reset();
    }

    constexpr void push(out_t sample){
        _history[_curr] = sample;
        _curr = (_curr + 1) % _arm_size;
    }
    constexpr out_t operate(size_t idx){
        out_t accum{};
        size_t index = _curr;
        for(size_t i = 0; i < _arm_size; ++i){
            accum += _taps[idx][i] * _history[index];
            index = index == 0 ? _arm_size - 1 : index - 1;
        }
        return accum;
    }
    constexpr out_t* operate_n(size_t idx, out_t* input, size_t n_in){
        out_t* out = new out_t[n_in];
        for(size_t i = 0; i < n_in; ++i){
            push(input[n_in]);
            out[i] = operate(idx);
        }
        return out;
    }

    constexpr void reset(){
        _curr = 0;
        std::fill_n(_history, _arm_size, out_t{});
    }
};

template<typename base_t, typename out_t>
class fir_resampler{
private:
    size_t P;
    size_t Q;
    base_t _scale;

    fir_pfb<base_t, out_t> _bank;
public:
    constexpr fir_resampler(size_t interp, size_t decim, base_t* taps, size_t num_taps){
        P = interp;
        Q = decim;
        _scale = P < Q ? std::sqrt(static_cast<base_t>(P) / static_cast<base_t>(Q)) : 1.0 / std::sqrt(static_cast<base_t>(P) / static_cast<base_t>(Q));
        _bank.update(P, taps, num_taps);
    }

    constexpr fir_pfb<base_t, out_t>* bank() noexcept{return &_bank;}

    constexpr void operate(out_t* samples, out_t* output){
        // Runs the resampler on a sample array of size Qx1
        // Output placed into array of size Px1
        size_t idx = 0;
        size_t n = 0;
        for(size_t i = 0; i < Q; ++i){
            _bank.push(samples[i]);

            while(idx < P){
                output[n++] = _bank.operate(idx) * _scale;
                idx += Q;
            }
            idx -= P;
        }
    }
    constexpr out_t* operate_n(out_t* samples, size_t n_in, size_t& n_written){
        // need to find the closest multiple for Q of the input size
        size_t pad_size = (Q - (n_in % Q)) % Q;
        size_t output_size = pad_size == 0 ? P * n_in / Q : std::floor(n_in * pad_size * P / Q);
        size_t iters = (n_in + pad_size) / Q;
        out_t* out = new out_t[output_size];
        out_t* temp;
        if(pad_size != 0){
            temp = new out_t[n_in + pad_size];
            std::uninitialized_default_construct_n(temp, n_in + pad_size);
            std::copy_n(samples, n_in, temp);
        } 
        else{
            temp = samples;
        }
        for(size_t i = 0; i < iters; ++i){
            operate(temp + i * Q, out + i * P);
        }
        n_written = output_size;
        return out;
    }

    constexpr void update(size_t interp, size_t decim, base_t* taps, size_t num_taps){
        P = interp;
        Q = decim;
        _scale = P < Q ? std::sqrt(static_cast<base_t>(P) / static_cast<base_t>(Q)) : 1.0 / std::sqrt(static_cast<base_t>(P) / static_cast<base_t>(Q));
        _bank.update(P, taps, num_taps);
    }

    constexpr void reset(){
        _bank.reset();
    }
};

template<typename base_t>
requires std::is_floating_point_v<base_t>
class digital_modem{
private:
    using cpx_t = std::complex<base_t>;
    static constexpr double PI = std::numbers::pi_v<base_t>;
    size_t _sps;
    constellation* _constel;
    base_t _loop_bw;
    base_t _rolloff;

    constexpr base_t sinc(base_t x){
        if(x > -1e-6 && x < 1e-6) return 1.0;
        else return std::sin(PI * x) / (PI * x);
    }

    /* Modulation-Demodulation Resampling */
    base_t _scale;
    fir_pfb<base_t, cpx_t> _upsample;
    fir_pfb<base_t, cpx_t> _downsample;

    constexpr void update_resample(){
        // Create two FIR root-raised cosine filters for the modulation and demodulation.
        // a filter size of 2 * sps * 32 + 1 works well for some reason (????)
        auto prototype = root_nyquist(_sps, _sps, 1.0, _rolloff, 2 * _sps * 32 + 1);
        _upsample.update(_sps, prototype.data(), prototype.size());
        _downsample.update(1, prototype.data(), prototype.size());
        _scale = 1.0 / std::sqrt(static_cast<base_t>(_sps));
    }
    /////////////////////////////////////////////
    
    /* Automatic Gain Control */
    base_t _agc_rate;
    base_t _agc_ref;
    base_t _agc_gain;
    base_t _agc_max_gain;

    constexpr void update_agc(){
        _agc_rate = 1e-4;
        _agc_ref = 1.0;
        _agc_gain = 1.0;
        _agc_max_gain = 1e3;
    }
    /////////////////////////////////////////////

    /* Preamble Energy Detector */
    uint8_t _preamble[8] = {0xAA, 0x0B, 0xF9, 0x42, 0x01, 0x2A, 0x18, 0xAF}; 
    size_t _ped_preamble_size;
    size_t* _ped_preamble;

    cpx_t* _ped_sync_words;
    size_t _ped_sync_size;
    fir<cpx_t, cpx_t> _ped_filter;

    base_t _ped_threshold;
    size_t _ped_corr_est;
    base_t _ped_time_est, _ped_phase_est;
    size_t _ped_expected_delay;

    base_t* _ped_hist;
    size_t _ped_curr;

    constexpr void update_ped(){
        size_t k = 8 / _constel->get_bps();
        _ped_preamble_size = 8 * k;
        _ped_preamble = new size_t[_ped_preamble_size];
        for(size_t i = 0; i < 8; ++i){
            uint8_to_point(_preamble[i], _ped_preamble + (i * k));
        }

        _ped_sync_size = _ped_preamble_size * _sps;
        _ped_sync_words = new cpx_t[_ped_sync_size];
        size_t m; // actual words written, might be different, not sure why though
        upsample_n(_ped_preamble, _ped_preamble_size, _ped_sync_words, m);
        for(size_t i = 0; i < _ped_sync_size; ++i){
            _ped_sync_words[i] = std::conj(_ped_sync_words[i]);
        }
        std::reverse(_ped_sync_words, _ped_sync_words + _ped_sync_size);
        _ped_filter.set_taps(_ped_sync_words, _ped_sync_size);

        base_t accum = 0.0;
        for(size_t i = 0; i < _ped_sync_size; ++i){
            accum += std::abs(_ped_sync_words[i] * std::conj(_ped_sync_words[i]));
        }
        _ped_threshold = 0.9 * std::pow(accum, 2.0);
        _ped_expected_delay = _ped_sync_size; // ideal point should be at the end of the preamble

        _ped_corr_est = 0;
        _ped_time_est = 0.0;
        _ped_phase_est = 0.0;

        _ped_hist = new base_t[3];
        _ped_curr = 0;
    }
    /////////////////////////////////////////////

    /* Phase Locked Loop */
    base_t _pll_phase, _pll_freq;
    base_t _pll_max_freq, _pll_min_freq;
    base_t _pll_damping, _pll_alpha, _pll_beta;
    base_t _pll_err;

    constexpr void update_pll(){
        _pll_phase = 0.0;
        _pll_freq = 0.0;
        _pll_max_freq = 0.25;
        _pll_min_freq = -0.25;
        _pll_damping = std::sqrt(2.0) / 2.0;
        _pll_alpha = 4.0 * _pll_damping * _loop_bw / (1.0 + 2.0 * _pll_damping * _loop_bw + std::pow(_loop_bw, 2.0));
        _pll_beta = 4.0 * std::pow(_loop_bw, 2.0) / (1.0 + 2.0 * _pll_damping * _loop_bw + std::pow(_loop_bw, 2.0));
        _pll_err = 0.0;
    }
    /////////////////////////////////////////////

    /* Frequency Locked Loop */
    base_t _fll_phase, _fll_freq;
    base_t _fll_max_freq, _fll_min_freq;
    base_t _fll_damping, _fll_alpha, _fll_beta;
    base_t _fll_err;
    cpx_t _fll_low, _fll_high;
    size_t _fll_size;
    fir<cpx_t, cpx_t> _fll_lower, _fll_upper;

    constexpr void update_fll(){
        _fll_damping = std::sqrt(2.0) / 2.0;
        _fll_alpha = 0.0;
        _fll_beta = 4.0 * _loop_bw / static_cast<base_t>(_sps);
        _fll_max_freq = 2.0 * PI / static_cast<base_t>(_sps);
        _fll_min_freq = -2.0 * PI / static_cast<base_t>(_sps);
        _fll_phase = 0.0;
        _fll_freq = 0.0;
        _fll_err = 0.0;
        _fll_low = 0.0;
        _fll_high = 0.0;
        _fll_size = 64;
        
        const base_t f_sps = static_cast<base_t>(_sps);
        size_t M = std::rint(static_cast<base_t>(_fll_size) / f_sps);
        base_t power = 0.0;
        base_t* temp = new base_t[_fll_size];
        for(size_t i = 0; i < _fll_size; ++i){
            base_t k = (static_cast<base_t>(i) * 2.0 / f_sps) - static_cast<base_t>(M);
            base_t pos = _rolloff * k;
            base_t tap = sinc(pos + 0.5) + sinc(pos - 0.5);
            power += std::pow(tap, 2.0);
            temp[i] = tap;
        }

        cpx_t* upper = new cpx_t[_fll_size];
        cpx_t* lower = new cpx_t[_fll_size];
        size_t N = (_fll_size - 1) / 2;
        for(size_t i = 0; i < _fll_size; ++i){
            base_t tap = temp[i] / power;
            base_t k = (static_cast<base_t>(i) - static_cast<base_t>(N)) * 0.5 / f_sps;
            size_t idx = _fll_size - i - 1;
            lower[idx] = tap * std::polar(1.0, -2.0 * PI * (1.0 + _rolloff) * k);
            upper[idx] = std::conj(lower[idx]);
        }
        _fll_lower.set_taps(lower, _fll_size);
        _fll_upper.set_taps(upper, _fll_size);
    }
    /////////////////////////////////////////////
public:
    // Convert a uint8 to a constellation point
    // points = size_t* size of [8 / bps]
    constexpr void uint8_to_point(uint8_t bit, size_t* points){
        size_t bps = _constel->get_bps();
        size_t k = 8 / bps;
        uint8_t mask = 0xFF >> (8 - bps);
        for(size_t i = 0; i < k; ++i){
            uint8_t new_mask = mask << (i * bps);
            uint8_t x = (bit & new_mask) >> (i * bps);
            points[i] = static_cast<size_t>(x);
        }
    }
    // Convert an array of points size 8 / bps to a single uint8
    constexpr void point_to_uint8(size_t* point, uint8_t& bit){
        size_t bps = _constel->get_bps();
        size_t k = 8 / bps; // expected array size
        bit = 0x00;
        for(size_t i = 0; i < k; ++i){
            std::bitset<8> U(point[i]); // U contains the bits 
            bit |= (U << (i * bps)).to_ullong();
        }
    }
public:
    digital_modem() = delete;
    constexpr digital_modem(const size_t& sps, const base_t& loop_bandwidth, const base_t& rolloff, constellation* const constel): _sps(sps), _constel(constel), _loop_bw(loop_bandwidth), _rolloff(rolloff){
        update();
    }

    constexpr void upsample(size_t point, cpx_t* output){
        cpx_t sample = _constel->get_point(point % _constel->get_size());
        size_t idx = 0; size_t n = 0;
        _upsample.push(sample);
        while(idx < _sps){
            output[n++] = _upsample.operate(idx) * _scale;
            idx += 1;
        }
        idx -= _sps;
    }
    constexpr void upsample_n(size_t* points, size_t n_in, cpx_t* output, size_t& n_written){
        size_t pad_size = (1 - (n_in % 1)) % 1;
        n_written = (n_in + pad_size) * _sps;
        size_t iters = n_in + pad_size;
        for(size_t i = 0; i < iters; ++i){
            if(i == iters - 1 && pad_size != 0) upsample(0, output + i * _sps);
            else upsample(points[i], output + i * _sps);
        }
    }

    constexpr void downsample(cpx_t* sample, size_t* output){
        size_t idx = 0; size_t n = 0;
        for(size_t i = 0; i < _sps; ++i){
            _downsample.push(sample[i]);
            while(idx < 1){
                output[n++] = _constel->decision(_downsample.operate(idx) * _scale);
                idx += _sps;
            }
            idx -= 1;
        }
    }
    constexpr void downsample_n(cpx_t* samples, size_t n_in, size_t* output, size_t& n_written){
        size_t pad_size = (_sps - (n_in % _sps)) % _sps;
        n_written = (n_in + pad_size) / _sps;
        size_t iters = (n_in + pad_size) / _sps;
        cpx_t* temp;
        if(pad_size != 0){
            temp = new cpx_t[n_in + pad_size];
            std::uninitialized_default_construct_n(temp, n_in + pad_size);
            std::copy_n(samples, n_in, temp);
        }
        else{
            temp = samples;
        }
        for(size_t i = 0; i < iters; ++i){
            downsample(temp + i * _sps, output + i);
        }
    }

    constexpr void downsample_cpx(cpx_t* sample, cpx_t* output){
        size_t idx = 0; size_t n = 0;
        for(size_t i = 0; i < _sps; ++i){
            _downsample.push(sample[i]);
            while(idx < 1){
                output[n++] = _downsample.operate(idx) * _scale;
                idx += _sps;
            }
            idx -= 1;
        }
    }
    constexpr void downsample_cpx_n(cpx_t* samples, size_t n_in, cpx_t* output, size_t& n_written){
        size_t pad_size = (_sps - (n_in % _sps)) % _sps;
        n_written = (n_in + pad_size) / _sps;
        size_t iters = (n_in + pad_size) / _sps;
        cpx_t* temp;
        if(pad_size != 0){
            temp = new cpx_t[n_in + pad_size];
            std::uninitialized_default_construct_n(temp, n_in + pad_size);
            std::copy_n(samples, n_in, temp);
        }
        else{
            temp = samples;
        }
        for(size_t i = 0; i < iters; ++i){
            downsample_cpx(temp + i * _sps, output + i);
        }
    }

    constexpr void agc(cpx_t sample, cpx_t* output){
        output[0] = sample * _agc_gain;
        _agc_gain += _agc_rate * (_agc_ref - std::norm(output[0]));
        _agc_gain = _agc_gain > _agc_max_gain ? _agc_max_gain : _agc_gain;
    }
    constexpr void agc_n(cpx_t* samples, size_t n_in, cpx_t* output, size_t& n_written){
        n_written = n_in;
        for(size_t i = 0; i < n_in; ++i){
            agc(samples[i], output + i);
        }
    }

    constexpr void fll(cpx_t sample, cpx_t* output){
        cpx_t nco = std::polar(1.0, _fll_phase);
        cpx_t fll_out = sample * nco;

        _fll_low = _fll_lower.filter_norm(fll_out);
        _fll_high = _fll_upper.filter_norm(fll_out);

        _fll_err = std::norm(_fll_high) - std::norm(_fll_low);
        _fll_freq += _fll_beta * _fll_err;
        _fll_phase += _fll_freq + _fll_alpha * _fll_err;

        if(_fll_phase > 2.0 * PI) _fll_phase = std::fmod(_fll_phase, 2.0 * PI);
        if(_fll_phase < -2.0 * PI) _fll_phase = std::fmod(_fll_phase, -2.0 * PI);
        _fll_freq = _fll_freq > _fll_max_freq ? _fll_max_freq : _fll_freq;
        _fll_freq = _fll_freq < _fll_min_freq ? _fll_min_freq : _fll_freq;

        output[0] = fll_out;
    }
    constexpr void fll_n(cpx_t* samples, size_t n_in, cpx_t* output, size_t& n_written){
        for(size_t i = 0; i < n_in; ++i){
            cpx_t nco = std::polar(1.0, _fll_phase);
            output[i] = samples[i] * nco;

            // _fll_lower.filter(output[i], &_fll_low);
            // _fll_upper.filter(output[i], &_fll_high);
            _fll_low = _fll_lower.filter_norm(output[i]);
            _fll_high = _fll_upper.filter_norm(output[i]);

            _fll_err = std::norm(_fll_high) - std::norm(_fll_low);
            _fll_freq += _fll_beta * _fll_err;
            _fll_phase += _fll_freq + _fll_alpha * _fll_err;

            if(_fll_phase > 2.0 * PI) _fll_phase = std::fmod(_fll_phase, 2.0 * PI);
            if(_fll_phase < -2.0 * PI) _fll_phase = std::fmod(_fll_phase, -2.0 * PI);
            _fll_freq = _fll_freq > _fll_max_freq ? _fll_max_freq : _fll_freq;
            _fll_freq = _fll_freq < _fll_min_freq ? _fll_min_freq : _fll_freq;
        }
        n_written = n_in;
    }

    constexpr void append_n(uint8_t* bits, size_t n_in, size_t* output, size_t& n_written){
        // size_t n = 0;
        size_t k = 8 / _constel->get_bps(); // number of points per byte
        size_t* points_buffer = new size_t[k];
        size_t n = 0;
        for(size_t i = 0; i < _ped_preamble_size; ++i){
            // this first batch is just our preamble points
            // TEST APPEND: IGNORE THIS SECTION
            output[n++] = _ped_preamble[i];  
        }
        for(size_t i = 0; i < n_in; ++i){
            // this second batch is our bits that need to be converted to a point
            uint8_to_point(bits[i], points_buffer);
            for(size_t j = 0; j < k; ++j){
                output[n++] = points_buffer[j];
            }
        }
        n_written = n;
    }

    constexpr void detect_n(cpx_t* samples, size_t n_in, cpx_t* output, size_t& n_written){
        size_t n = 0;
        bool detected = false;
        bool write = false;
        for(size_t i = 0; i < n_in; ++i){
            _ped_filter.insert(samples[i]);
            cpx_t ped_corr = _ped_filter.operate();
            _ped_filter.increment();

            base_t corr_mag = std::norm(ped_corr);
            _ped_hist[_ped_curr] = corr_mag;
            _ped_curr = (_ped_curr + 1) % 3;
            if(corr_mag >= _ped_threshold){
                // we have a good correlation, so begin writing to the output, and make some estimates
                detected = true;
                write = true;
            }
            
            if(detected){
                // on the first pass once detected, make estimates, and then continue writing (skip the preamble in output data)
                _ped_corr_est = i; // point where correlation happens
                _ped_phase_est = std::atan2(ped_corr.imag(), ped_corr.real());
                size_t prev = _ped_curr == 0 ? 2 : _ped_curr - 1;
                size_t next = (_ped_curr + 1) % 3;
                base_t numer = _ped_hist[prev] + 2.0 * _ped_hist[_ped_curr] + 3.0 * _ped_hist[next];
                base_t denom = _ped_hist[prev] + _ped_hist[_ped_curr] + _ped_hist[next];
                _ped_time_est = numer / denom - 2.0;
                detected = false;
                continue;
            }
            if(write){
                output[n++] = samples[i];
            }
        }
        n_written = n;
    }

    constexpr void pll_n(cpx_t* input, size_t n_in, cpx_t* output, size_t& n_written){
        bool ped = true; // we will use the ped detection initially, and then use our own phase detector
        for(size_t i = 0; i < n_in; ++i){
            cpx_t nco;
            if(ped){
                nco = std::polar(1.0, -_ped_phase_est);
                ped = false;
            }
            else{
                nco = std::polar(1.0, -_pll_phase);
            }
            output[i] = input[i] * nco;

            _pll_err = _constel->phase_error_detector(output[i]);

            _pll_freq += _pll_beta * _pll_err;
            _pll_phase += _pll_freq + _pll_alpha * _pll_err;

            if(_pll_phase > 2.0 * PI) _pll_phase = std::fmod(_pll_phase, 2.0 * PI);
            if(_pll_phase < -2.0 * PI) _pll_phase = std::fmod(_pll_phase, -2.0 * PI);
            _pll_freq = _pll_freq > _pll_max_freq ? _pll_max_freq : _pll_freq;
            _pll_freq = _pll_freq < _pll_min_freq ? _pll_min_freq : _pll_freq;
        }
        n_written = n_in;
    }

    /**
     * @brief Primary modulation function. Performs conversion from bit stream to complex symbols.
     * @param bits input bit stream
     * @param n_in number of bits in stream
     * @param output complex output symbols
     * @param n_written number of symbols written to the output
     */
    constexpr void modulate_n(uint8_t* bits, size_t n_in, cpx_t* output, size_t& n_written){
        // loop over n_in + _ped_preamble_size (total number of points in the message)
        // for(size_t i = 0; )
        size_t k = 8 / _constel->get_bps(); // number of points per byte
        size_t* points_buffer = new size_t[k];
        size_t n = 0;
        for(size_t i = 0; i < _ped_preamble_size; ++i){
            // this first batch is just our preamble points
            // output[n++] = _ped_preamble[i];  
            upsample(_ped_preamble[i], output + i * _sps);
            n += _sps;
        }
        // n += _sps;
        // we got _ped_preamble_size * sps samples last iteration
        for(size_t i = 0; i < n_in; ++i){
            // this second batch is our bits that need to be converted to a point
            uint8_to_point(bits[i], points_buffer); // data held in the points buffer
            for(size_t j = 0; j < k; ++j){
                // output[n++] = points_buffer[j];
                upsample(points_buffer[j], output + n);
                n += _sps;
            }
        }
        n_written = n;
    }

    /**
     * @brief Demodulates a complex stream to a stream of downsampled complex symbols
     * @param input input complex symbols
     * @param n_in number of complex symbols
     * @param output bit stream output
     * @param n_written number of bits written to the output
     */
    constexpr void demodulate_cpx_n(cpx_t* input, size_t n_in, cpx_t* output, size_t& n_written){
        cpx_t* buffer = new cpx_t[n_in]; // main storage buffer
        // agc_n(input, n_in, buffer, n_written);
        // fll_n(buffer, n_written, buffer, n_written);
        fll_n(input, n_in, buffer, n_written);
        detect_n(buffer, n_written, buffer, n_written);
        pll_n(buffer, n_written, buffer, n_written);
        downsample_cpx_n(buffer, n_written, output, n_written);
    }

    /**
     * @brief Primary demodulation function. Performs syncing and equalization from complex symbols to a bit stream.
     * @param input 
     * @param n_in 
     * @param bits 
     * @param n_written 
     */
    constexpr void demodulate_n(cpx_t* input, size_t n_in, uint8_t* bits, size_t& n_written){
        size_t bps = _constel->get_bps();
        size_t k = 8 / bps;
        n_written = 0;
        
        cpx_t* cpx_buffer = new cpx_t[n_in]; // internal storage buffer for holding complex samples

        // Detection, Equalization
        fll_n(input, n_in, cpx_buffer, n_written); // rough carrier centering
        detect_n(cpx_buffer, n_written, cpx_buffer, n_written); // preamble detector, removes the preamble entirely
        pll_n(cpx_buffer, n_written, cpx_buffer, n_written); // fine phase estimation of symbols

        // Downsampling to Constellation
        size_t* points_buffer = new size_t[n_written / _sps]; // internal storage buffer for holding constellation points
        downsample_n(cpx_buffer, n_written, points_buffer, n_written); // downsamples the recovered data into constellation points

        std::cout<<"Downsample Wrote: "<<n_written<<"\n";

        // Reconvert to a bitstream
        size_t* temp = new size_t[k]; // must hold k points to make a single uint8
        size_t iters = n_written / k; // number of input data bits (probably)
        std::cout<<"iters: "<<iters<<"\n";
        for(size_t i = 0; i < iters; ++i){
            // load into the temp
            for(size_t j = 0; j < 8; ++j){
                temp[j] = points_buffer[i * k + j];
            }
            point_to_uint8(temp, bits[i]);
        }
        n_written = iters;
    }

    /**
     * @brief Decode an array of points into the binary data it represented as uint8's.
     * @param input 
     * @param n_in 
     * @param output 
     * @param n_written 
     */
    constexpr void decode_n(size_t* input, size_t n_in, uint8_t* output, size_t& n_written){
        size_t bps = _constel->get_bps();
        size_t k = 8 / bps; // number of points we take per uint8
        size_t* points_buffer = new size_t[k];
        std::uninitialized_default_construct_n(points_buffer, k); 
    }

    constexpr size_t sps() const noexcept{return _sps;}
    constexpr constellation* constel() const noexcept{return _constel;}
    constexpr base_t loop_bandwidth() const noexcept{return _loop_bw;}
    constexpr fir_pfb<base_t, cpx_t>* upsampler() const noexcept{return &_upsample;}
    constexpr fir_pfb<base_t, cpx_t>* downsampler() const noexcept{return &_downsample;}
    constexpr size_t preamble_len() const noexcept{return _ped_preamble_size;}
    constexpr size_t* preamble() const noexcept{return _ped_preamble;}
    constexpr size_t correlation_est() const noexcept{return _ped_corr_est;}
    constexpr base_t time_est() const noexcept{return _ped_time_est;}
    constexpr base_t phase_est() const noexcept{return _ped_phase_est;}
    constexpr fir<cpx_t, cpx_t> fll_lower() const noexcept{return _fll_lower;}
    constexpr fir<cpx_t, cpx_t> fll_upper() const noexcept{return _fll_upper;}

    constexpr void set_sps(const size_t& sps){_sps = sps; update();}
    constexpr void set_constellation(constellation* const constel){_constel = constel; update();}
    constexpr void set_loop_bandwidth(const base_t& loop_bandwidth){_loop_bw = loop_bandwidth; update();}

    constexpr void update(){
        update_resample();
        update_agc();
        update_ped();
        update_pll();
        update_fll();
        reset();
    }
    constexpr void reset(){
        _upsample.reset();
        _downsample.reset();
        _ped_filter.reset();
        _fll_lower.reset();
        _fll_upper.reset();
    }
};

#endif