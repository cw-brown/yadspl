#ifndef TX_H
#define TX_H

#include <vector>
#include <complex>
#include <numbers>
#include <numeric>
#include <cmath>

#include "polyphase.hpp"
#include "constellations.hpp"

/**
 * @brief Constellation Modulator is a transmitter that takes in a bit stream and returns complex values
 * 
 */
class constellation_modulator{
private:
    size_t _sps;
    double _alpha;
    constellation _constel;

    std::vector<double> _taps;
public:
    constellation_modulator(constellation constellation, const size_t& sps){}


};

#endif