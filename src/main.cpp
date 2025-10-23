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
#include "modem.hpp"

#include "fft.hpp"

#define DO_WINDOW true

int main(){
    using cpx = std::complex<double>;
    size_t sps = 4;
    size_t n_filt = 32;
    size_t points = 250;
    size_t N = sps * points;

    static float v = 0.0f;
    static float offset = 0.0f;

    constellation_16qam constel{};
    noise<double> sig_gen{};
    channel_model channel(offset, v);
    rectangular_modulator modulator(&constel, sps, n_filt, 0.35);
    firefighter recovery(&constel, sps, n_filt, 2.0*3.1415/100.0, 0.35);

    cpx* data = new cpx[N];
    cpx* data_copy = new cpx[N];
    for(size_t i = 0; i < points; ++i){
        auto point = sig_gen.random_int_range(0, constel.get_size() - 1);
        modulator.operate(point, data + i * sps);
    }
    double* data_fft = new double[N];
    double* data_freq = new double[N];

    cpx* output = new cpx[N];
    double* output_fft = new double[N];
    double* output_freq = new double[N];

    std::ring<double> real(points);
    std::ring<double> imag(points);
    cpx* const_buff = new cpx;
    auto prot = root_nyquist(n_filt, n_filt * sps, 1.0, 0.35, 8 * sps * n_filt);
    resampler<cpx> arb(1.0 / sps, n_filt, prot);

#if DO_WINDOW
    GLFWwindow* window = glfw_makeNewWindow(1920, 1080, "Yet Another DSP Library", true, true, true);
    ImPlot::CreateContext();
    glfwSetKeyCallback(window, key_call);

    while(!glfwWindowShouldClose(window)){
        glfw_frame();

        std::copy_n(data, N, data_copy);
        for(size_t i = 0; i < N; ++i){
            channel.operate(data_copy + i);
        }
        compute_psd(data_copy, N, data_fft, data_freq);
        for(size_t i = 0; i < N; ++i){
            recovery.operate(data_copy[i], output + (int)i);
        }
        compute_psd(output, N, output_fft, output_freq);
        for(size_t i = 0; i < N; ++i){
            if(arb.operate(output[i], const_buff) != 0){
                real.push_back(const_buff->real());
                imag.push_back(const_buff->imag());
            }
            
        }


        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x, ImGui::GetIO().DisplaySize.y));
        ImGuiWindowFlags topbarflags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize 
            | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoCollapse;
    
        ImGui::Begin("Plottings", nullptr, topbarflags);
        if(ImGui::BeginTabBar("Main Tabs")){
        if(ImGui::BeginTabItem("Item 0")){
            ImGui::SliderFloat("Offset", &offset, -1.0, 1.0);
            ImGui::SliderFloat("Noise", &v, 0.0, 1.0);
            channel.set_noise(v);
            channel.set_offset(offset);
            ImGui::Text((std::string("Lower: ") + std::to_string(recovery.fll_upper())).c_str());
            ImGui::Text((std::string("Upper: ") + std::to_string(recovery.fll_lower())).c_str());
            if(ImPlot::BeginPlot("Plot 0", ImVec2(-1, 750))){
                ImPlot::SetupAxisLimits(ImAxis_Y1, -30, 5);
                ImPlot::PlotLine("Data", data_freq, data_fft, N);
                ImPlot::PlotLine("Output", output_freq, output_fft, N);
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginSubplots("Debug Info", 1, 2, ImVec2(-1, 750))){
            if(ImPlot::BeginPlot("Errors")){
                ImPlot::SetupAxesLimits(0, 500, 2.0*3.141, -2.0*3.141);
                ImPlot::PlotLine("Error", recovery.fll_err(), 500);
                ImPlot::PlotLine("Phase", recovery.fll_phase(), 500);
                ImPlot::PlotLine("Frequency", recovery.fll_freq(), 500);
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Constellation")){
                ImPlot::SetupAxesLimits(-2, 2, -2, 2);
                ImPlot::PlotScatter("Constellation", real.data(), imag.data(), real.size());
                ImPlot::EndPlot();
            }
                ImPlot::EndSubplots();
            }

            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("Items 1")){

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

