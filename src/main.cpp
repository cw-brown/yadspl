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

#include <fftw3.h>

#include "implot.h"

#include "helper_funcs.h"
#include "fir_filter.hpp"
#include "noise.hpp"
#include "yadpsl_math.hpp"
#include "constellations.hpp"
#include "symbol_rec.hpp"
#include "ring.hpp"
#include "polyphase.hpp"
#include "graphics.hpp"
#include "modem2.hpp"

#include "fft.hpp"

#define DO_WINDOW false


int main(){
    size_t sps = 8;
    static size_t n = 111;
    static float nois = 0.0;
    static float offset = 0.0;
    static noise<double> sig_gen{};
    static constellation_bpsk constel{};
    static channel_model channel(offset, nois);
    digital_modem<double> modem(sps, 2.0*3.1/100.0, 0.35, &constel);
    size_t m = modem.preamble_len();

    uint8_t* test_data = new uint8_t[n]; 
    for(size_t i = 0; i < n; ++i){
        test_data[i] = static_cast<uint8_t>(sig_gen.random_int_range(0, 255));
    }

    size_t total = m + n * (8 / constel.get_bps());
    size_t* points = new size_t[total]; // time data with preamble
    size_t k;
    modem.append_n(test_data, n, points, k);
    size_t* indx = new size_t[k];
    std::iota(indx, indx + k, 0);

    size_t expected = (m * sps) + (n * 8 / constel.get_bps()) * sps; // expected number, actual might be lower
    // all buffers can be that size, since nothing increases the total size, only decreases
    std::complex<double>* tx_data_buffer = new std::complex<double>[expected];
    std::complex<double>* rx_data_buffer = new std::complex<double>[expected]; 
    size_t* data_buffer = new size_t[expected];
    uint8_t* output_data = new uint8_t[n + 300];
    std::uninitialized_default_construct_n(tx_data_buffer, expected);
    std::uninitialized_default_construct_n(rx_data_buffer, expected);
    
    size_t tx_n = 0;
    size_t rx_n = 0;
    size_t bits_n = 0;

    double* input_fft = new double[expected];
    double* input_freq = new double[expected];
    double* output_fft = new double[expected];
    double* output_freq = new double[expected];
    double* real_constellation = new double[expected];
    double* imag_constellation = new double[expected];

    std::cout<<"Premable Len: "<<m<<", N Data: "<<n<<", Total Const Points: "<<k<<"\n";
    std::cout<<"Expected Upsample Size: "<<expected<<"\n";

    // size_t w[8];
    // modem.uint8_to_point(0xFE, w);
    // for(size_t i = 0; i < 8; ++i){
    //     std::cout<<"Point "<<i<<": "<<w[i]<<"\n";
    // }

    // uint8_t back = 0x00;
    // modem.point_to_uint8(w, back);
    // std::cout<<std::hex<<"Back Converted: "<<static_cast<int>(back)<<"\n";

    modem.modulate_n(test_data, n, tx_data_buffer, tx_n);

    modem.demodulate_n(tx_data_buffer, tx_n, output_data, rx_n);
    modem.demodulate_n(tx_data_buffer, tx_n, output_data, rx_n);
    modem.demodulate_n(tx_data_buffer, tx_n, output_data, rx_n);
    modem.demodulate_n(tx_data_buffer, tx_n, output_data, rx_n);
    modem.demodulate_n(tx_data_buffer, tx_n, output_data, rx_n);
    std::cout<<"RX Wrote: "<<rx_n<<"\n";

#if DO_WINDOW
    GLFWwindow* window = glfw_makeNewWindow(1920, 1080, "Yet Another DSP Library", true, true, true);
    ImPlot::CreateContext();
    glfwSetKeyCallback(window, key_call);

    while(!glfwWindowShouldClose(window)){
        glfw_frame();

        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x, ImGui::GetIO().DisplaySize.y));
        ImGuiWindowFlags topbarflags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize 
            | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoCollapse;
    
        ImGui::Begin("Plottings", nullptr, topbarflags);
        if(ImGui::BeginTabBar("Main Tabs")){
        if(ImGui::BeginTabItem("Item 0")){
            static int dat_off = 0;
            ImGui::SliderFloat("Noise", &nois, 0.0, 0.5);
            ImGui::SliderFloat("Offset", &offset, 0.0, 1.0);
            ImGui::SliderInt("Data Offset", &dat_off, 0, 64);
            channel.set_noise(nois);
            channel.set_offset(offset);
            std::string display = std::string("Time Estimate: ") + std::to_string(modem.time_est()) + 
                                std::string(", Phase Estimate: ") + std::to_string(modem.phase_est()) + 
                                std::string(", Point Estimate: ") + std::to_string(modem.correlation_est());
                                
            ImGui::Text(display.c_str());

            modem.modulate_n(test_data, n, tx_data_buffer, tx_n);
            for(size_t i = 0; i < tx_n; ++i){
                channel.operate(tx_data_buffer + i); // add freq offset and noise
            }
            compute_psd(tx_data_buffer, tx_n, input_fft, input_freq);

            // modem.fll_n(tx_data_buffer, tx_n, rx_data_buffer, rx_n);
            // size_t r = rx_n;
            // compute_psd(rx_data_buffer, r, output_fft, output_freq);

            // modem.detect_n(rx_data_buffer, rx_n, rx_data_buffer, rx_n);
            // std::string pedstuff = std::string("PED Wrote: ") + std::to_string(rx_n);
            // ImGui::Text(pedstuff.c_str());

            // modem.pll_n(rx_data_buffer, rx_n, rx_data_buffer, rx_n);
            // modem.downsample_cpx_n(rx_data_buffer, rx_n, rx_data_buffer, rx_n);

            // // modem.demodulate_n(tx_data_buffer, tx_n, rx_data_buffer, rx_n);
            // for(size_t i = 0; i < rx_n; ++i){
            //     real_constellation[i] = rx_data_buffer[i].real();
            //     imag_constellation[i] = rx_data_buffer[i].imag();
            //     data_buffer[i] = constel.decision(rx_data_buffer[i]);
            // }

            modem.demodulate_n(tx_data_buffer, tx_n, output_data, bits_n);

            if(ImPlot::BeginSubplots("Data", 1, 2, ImVec2(-1, 700))){
            if(ImPlot::BeginPlot("Spectrum")){
                ImPlot::PlotLine("Input", input_freq, input_fft, tx_n);
                // ImPlot::PlotLine("Output", output_freq, output_fft, r);
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Time Data")){
                // ImPlot::PlotStairs("Preamble", indx, points, m);
                // ImPlot::PlotStairs("Data", indx + m, points + m, k - m);
                // ImPlot::PlotStairs("Output", indx + dat_off, data_buffer, rx_n - dat_off);
                ImPlot::PlotStairs("Test Data", test_data, n);
                ImPlot::PlotStairs("Received", output_data, bits_n);
                ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
            }
            if(ImPlot::BeginSubplots("Debug", 1, 2, ImVec2(-1, 700))){
            if(ImPlot::BeginPlot("Constellation")){
                ImPlot::SetupAxesLimits(-2, 2, -2, 2);
                // ImPlot::PlotScatter("Constellation", real_constellation, imag_constellation, rx_n);
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Errors")){
                ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
            }
            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("Items 1")){
            if(ImPlot::BeginPlot("Plot 10", ImVec2(-1, 750))){
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Plot 1", ImVec2(-1, 750))){
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("Items 2")){
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
        }
        ImGui::End();

        glfw_render(window);
    }

    glfw_cleanup(window);
    ImPlot::DestroyContext();
#endif
    return 0;
}

