#include <iostream>
#include <complex>
#include <array>
#include <memory>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <random>
#include <algorithm>
#include <utility>


#include "implot.h"

#include "helper_funcs.h"
#include "fir_filter.hpp"
// #include "iir_filter.hpp"
// #include "polynomial.hpp"
#include "noise.hpp"
#include "yadpsl_math.hpp"
#include "constellations.hpp"
#include "symbol_rec.hpp"
#include "ring.hpp"
#include "transmitter.hpp"
#include "polyphase.hpp"
#include "graphics.hpp"
// #include "filter.hpp"

#include "fft.hpp"

void key_call(GLFWwindow* window, int key, int, int action, int){
    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
}

#define DO_WINDOW true

int main(){
    size_t sps = 4;
    double filter_bw = 0.35;
    size_t n_filts = 32;
    double loop_bw = 2.0*3.1415/100.0;

    auto tx_prot = root_nyquist(n_filts, n_filts, 1.0, filter_bw, 8*n_filts*sps);
    auto rx_prot = root_nyquist(n_filts, n_filts*sps, 1.0, filter_bw, 8*n_filts*sps);

    static constellation_qpsk constel{};
    static noise<double> sig_gen{};
    static polyphase_upsampler resampling_filter(sps, tx_prot, n_filts);

    // Setup a data array for constellation
    std::vector<std::complex<double>> constellation_data;
    std::vector<double> real_const;
    std::vector<double> imag_const;

    constellation_data.reserve(1024);
    real_const.reserve(1024);
    imag_const.reserve(1024);

    for(size_t i = 0; i < 1024; ++i){
        std::complex<double> point = constel.get_point(sig_gen.randomValue(constel.get_bps()));
        constellation_data.push_back(point);
        real_const.push_back(point.real() + 0.01*sig_gen.randomDouble());
        imag_const.push_back(point.imag() + 0.01*sig_gen.randomDouble());
    }

    // Setup a data array for the filtered symbols
    std::vector<std::complex<double>> tx_data = resampling_filter.filterN(constellation_data);
    std::vector<double> real_tx(tx_data.size());
    std::vector<double> imag_tx(tx_data.size());

    // Spectrum for tx data
    std::vector<std::complex<double>> tx_data_fft(tx_data);
    std::vector<double> tx_spectrum(tx_data.size(), 0.0);
    std::vector<double> tx_spectrum_freqs(tx_spectrum.size());
    for(size_t j = 0; j < tx_spectrum.size(); ++j){
        double half = std::round(tx_spectrum.size() / 2);
        tx_spectrum_freqs[j] = (static_cast<double>(j) - half)*1e-3;
    }

    // FLL edge recovery
    static frequency_recovery freq_rec(sps, filter_bw, loop_bw, n_filts);

    std::vector<std::complex<double>> rx_data = freq_rec.operate(tx_data);
    std::vector<std::complex<double>> rx_data_fft(rx_data);
    std::vector<double> rx_spectrum(rx_data.size(), 0.0);

    // PLL recovery
    static phase_recovery phase_rec(sps, loop_bw, n_filts, 1.5, filter_bw);

    std::vector<std::complex<double>> rec_data = phase_rec.operate(rx_data);
    std::vector<std::complex<double>> rec_data_fft(rec_data);
    std::vector<double> rec_spectrum(rec_data.size(), 0.0);
    std::vector<double> rx_spectrum_freqs(rec_spectrum.size(), 0.0);
    for(size_t i = 0; i < rx_spectrum_freqs.size(); ++i){
        double half = std::round(rx_spectrum_freqs.size() / 2);
        rx_spectrum_freqs[i] = i - half;
    }

    std::vector<double> real_rec(rec_data.size());
    std::vector<double> imag_rec(rec_data.size());
    for(size_t i = 0; i < rec_data.size(); ++i){
        real_rec[i] = rec_data[i].real();
        imag_rec[i] = rec_data[i].imag();
    }

    // Symbol recovery constellation
    static symbol_recovery symb_rec(&constel, loop_bw, -0.25, 0.25);
    std::vector<unsigned int> symb_const = symb_rec.operate(rec_data);
    std::vector<double> symb_real(symb_const.size());
    std::vector<double> symb_imag(symb_const.size());
    for(size_t i = 0; i < symb_const.size(); ++i){
        auto point = constel.get_point(symb_const[i]);
        symb_real[i] = point.real();
        symb_imag[i] = point.imag();
    }

    static bool overlay = false;
    size_t i = 0;

    if(DO_WINDOW){
    GLFWwindow* window = glfw_makeNewWindow(1920, 1080, "Yet Another DSP Library", true, true, true);
    ImPlot::CreateContext();
    glfwSetKeyCallback(window, key_call);

    while(!glfwWindowShouldClose(window)){
        glfw_frame();

        // Add more data for constellation
        std::complex<double> point = constel.get_point(sig_gen.randomValue(constel.get_bps()));
        constellation_data.erase(constellation_data.begin());
        real_const.erase(real_const.begin());
        imag_const.erase(imag_const.begin());
        const std::complex<double> awgn(0.01*sig_gen.randomDouble(), 0.01*sig_gen.randomDouble());
        constellation_data.push_back(point + awgn);
        real_const.push_back(point.real() + awgn.real());
        imag_const.push_back(point.imag() + awgn.imag());

        // Filter the constellation data
        tx_data = resampling_filter.filterN(constellation_data);
        for(size_t j = 0; j < tx_data.size(); ++j){
            real_tx[j] = tx_data[j].real();
            imag_tx[j] = tx_data[j].imag();
        }

        if(i % 64 == 0){
            tx_data_fft = tx_data;
            fft_transform(tx_data_fft, false);
            tx_spectrum = make_psd(tx_data_fft);

            rx_data = freq_rec.operate(tx_data);
            rx_data_fft = rx_data;
            fft_transform(rx_data_fft, false);
            rx_spectrum = make_psd(rx_data_fft);
        }

        if(i % 64 == 0){
                rec_data = phase_rec.operate(rx_data);
                rec_data_fft = rec_data;
                fft_transform(rec_data_fft, false);
                rec_spectrum = make_psd(rec_data_fft, false);

                symb_const = symb_rec.operate(rec_data);
        }
        i++;

        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x, ImGui::GetIO().DisplaySize.y));
        ImGuiWindowFlags topbarflags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize 
            | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoCollapse;
    
        ImGui::Begin("Plottings", nullptr, topbarflags);
        ImGui::BeginTabBar("Main Tabs");
        if(ImGui::BeginTabItem("TX Data")){
            if(ImGui::Button("Overlay Recovery", ImVec2(300, 50))) overlay = !overlay;
            if(ImPlot::BeginPlot("Spectrum", ImVec2(-1, 750))){
                ImPlot::SetupAxes("Frequency (kHz)", "Magnitude (dB)");
                ImPlot::SetupAxisLimits(ImAxis_Y1, -30, 5);
                ImPlot::PlotLine("TX Spectrum", tx_spectrum_freqs.data(), tx_spectrum.data(), tx_spectrum.size());
                if(overlay){
                    ImPlot::PlotLine("RX Spectrum", tx_spectrum_freqs.data(), rx_spectrum.data(), tx_spectrum.size());
                }
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginSubplots("", 1, 2, ImVec2(-1, 750))){
                if(ImPlot::BeginPlot("Constellation IQ Plot")){
                    ImPlot::SetupAxesLimits(-1.0, 1.0, -1.0, 1.0);
                    ImPlot::SetupAxes("In-Phase", "Quadrature");
                    ImPlot::PlotScatter("Constellation", real_const.data(), imag_const.data(), imag_const.size());
                    ImPlot::EndPlot();
                }
                if(ImPlot::BeginPlot("Transmitted Data")){
                    ImPlot::SetupAxes("Sample", "Amplitude");
                    ImPlot::SetupAxesLimits(64, 512, -1.45, 1.45);
                    ImPlot::SetNextLineStyle(ImVec4(0,0,0,-1), 4.0);
                    ImPlot::PlotLine("Real", real_tx.data(), 2048);
                    ImPlot::SetNextLineStyle(ImVec4(0,0,0,-1), 4.0);
                    ImPlot::PlotLine("Imaginary", imag_tx.data(), 2048);
                    ImPlot::EndPlot();
                }
                ImPlot::EndSubplots();
            }
            ImGui::EndTabItem();
        }
        
        if(ImGui::BeginTabItem("RX Data")){
            for(size_t j = 0; j < rec_data.size(); ++j){
                real_rec[j] = rec_data[j].real();
                imag_rec[j] = rec_data[j].imag();
            }
            for(size_t j = 0; j < symb_const.size(); ++j){
                auto point = constel.get_point(symb_const[j]);
                symb_real[j] = point.real();
                symb_imag[j] = point.imag();
            }
            if(ImPlot::BeginPlot("Recovered Spectrum", ImVec2(-1, 750))){
                ImPlot::SetupAxes("Frequency (kHz)", "Magnitude (dB)");
                ImPlot::SetupAxisLimits(ImAxis_Y1, -30, 5);
                ImPlot::PlotLine("", rx_spectrum_freqs.data(), rec_spectrum.data(), rec_spectrum.size());
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Recovered Constellation", ImVec2(-1, 750))){
                ImPlot::SetupAxesLimits(-1.0, 1.0, -1.0, 1.0);
                ImPlot::SetupAxes("In-Phase", "Quadrature");
                ImPlot::PlotScatter("", symb_real.data(), symb_imag.data(), symb_real.size());
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
        ImGui::End();

        glfw_render(window);
    }

    glfw_cleanup(window);
    ImPlot::DestroyContext();
    }
    return 0;
}

