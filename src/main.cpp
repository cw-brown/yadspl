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
    size_t sps = 2;
    size_t n_filt = 32;
    size_t points = 150;

    static float v = 0.0f;
    static float offset = 0.0f;

    constellation_bpsk constel{};
    noise<double> sig_gen{};
    channel_model channel(offset, v);
    rectangular_modulator modulator(&constel, sps, n_filt, 0.35);
    firefighter recovery(&constel, sps, n_filt, 2.0*3.1415/100.0, 0.35);

    size_t pre_len = recovery.get_preamble_size();
    size_t* preamble = recovery.get_preamble_points();
    size_t N = points + pre_len;
    size_t* bin_data = new size_t[N];
    size_t* bin_idx = new size_t[N];

    cpx* input = new cpx[N * sps];
    double* input_fft = new double[N * sps];
    double* input_freq = new double[N * sps];

    for(size_t i = 0; i < pre_len; ++i){
        bin_data[i] = preamble[i];
        modulator.operate(preamble[i], input + i * sps);
        bin_idx[i] = i;
    }
    for(size_t i = pre_len; i < N; ++i){
        size_t point = sig_gen.random_int_range(0, constel.get_size() - 1);
        bin_data[i] = point;
        modulator.operate(point, input + i * sps);
        bin_idx[i] = i;
    }
    compute_psd(input, N * sps, input_fft, input_freq);

    cpx* corr = new cpx[N * sps];
    double* corr_mag = new double[N * sps];
    for(size_t i = 0; i < N * sps; ++i){
        corr[i] = recovery.preamble(input[i]);
    }
    for(size_t i = 0; i < N * sps; ++i){
        corr_mag[i] = std::norm(corr[i]);
    }

    for(size_t i = 0; i < N * sps; ++i){
        auto b = recovery.pev(input[i]);
        if(b){break;}
    }

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
            if(ImPlot::BeginSubplots("Data", 1, 2, ImVec2(-1, 750))){
            if(ImPlot::BeginPlot("Spectrum")){
                ImPlot::PlotLine("", input_freq, input_fft, N * sps);
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Time Data")){
                ImPlot::PlotLine("Preamble", bin_idx, bin_data, pre_len);
                ImPlot::PlotLine("Data", bin_idx + pre_len, bin_data + pre_len, points);
                ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
            }
            if(ImPlot::BeginSubplots("Debug", 1, 2, ImVec2(-1, 750))){
            if(ImPlot::BeginPlot("Constellation")){
                ImPlot::PlotLine("", corr_mag, N * sps);
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
            if(ImPlot::BeginPlot("Plot 0", ImVec2(-1, 750))){
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

