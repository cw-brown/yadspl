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

#define DO_WINDOW true

int main(){
    size_t sps = 4;
    static size_t n = 250;
    static float nois = 0.0;
    static float offset = 0.0;
    static noise<double> sig_gen{};
    static constellation_bpsk constel{};
    static channel_model channel(offset, nois);
    digital_modem<double> modem(sps, 2.0*3.1/100.0, &constel);
    size_t m = modem.preamble_len();

    uint8_t* test_data = new uint8_t[n]; 
    // size_t* indices = new size_t[n]; std::iota(indices, indices + n, 0);
    for(size_t i = 0; i < n; ++i){
        test_data[i] = static_cast<uint8_t>(sig_gen.random_int_range(0, 255));
    }

    size_t total = m + n * (8 / constel.get_bps());
    size_t* points = new size_t[total];
    size_t k;
    modem.append_n(test_data, n, points, k);
    size_t* indx = new size_t[k];
    std::iota(indx, indx + k, 0);

    size_t expected = (m * sps) + (n * 8 / constel.get_bps()) * sps;
    std::complex<double>* buffer = new std::complex<double>[expected];
    double* fft = new double[expected];
    double* freq = new double[expected];

    size_t* outputs = new size_t[total];
    size_t o;

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
            ImGui::SliderFloat("Noise", &nois, 0.0, 0.5);
            ImGui::SliderFloat("Offset", &offset, 0.0, 1.0);
            channel.set_noise(nois);
            channel.set_offset(offset);
            ImGui::Text((std::string("Time Estimate: ") + std::to_string(modem.time_est())).c_str());
            ImGui::Text((std::string("Phase Estimate: ") + std::to_string(modem.phase_est())).c_str());

            size_t u;
            modem.modulate_n(test_data, n, buffer, u);
            for(size_t i = 0; i < expected; ++i){
                channel.operate(buffer + i);
            }
            compute_psd(buffer, u, fft, freq);

            modem.detect_n(buffer, expected, buffer, u);

            modem.downsample_n(buffer, expected, outputs, o);
    

            if(ImPlot::BeginSubplots("Data", 1, 2, ImVec2(-1, 750))){
            if(ImPlot::BeginPlot("Spectrum")){
                ImPlot::PlotLine("", freq, fft, u);
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Time Data")){
                ImPlot::PlotStairs("Preamble", indx, points, m);
                ImPlot::PlotStairs("Data", indx + m, points + m, k - m);
                ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
            }
            if(ImPlot::BeginSubplots("Debug", 1, 2, ImVec2(-1, 750))){
            if(ImPlot::BeginPlot("Constellation")){
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
            static int arm = 0;
            ImGui::SliderInt("Arm", &arm, 0, 1);
            if(ImPlot::BeginPlot("Plot 0", ImVec2(-1, 750))){
                // ImPlot::PlotBars("", bank->arm(arm), bank->arm_size());
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

