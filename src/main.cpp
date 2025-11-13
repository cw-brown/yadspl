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

    static size_t points = 1500;
    static noise<double> sig_gen{};
    static constellation_bpsk constel{};
    std::complex<double>* outs = new std::complex<double>[points];
    for(size_t i = 0; i < points; ++i){
        auto point = sig_gen.random_int_range(0, constel.get_size() - 1);
        outs[i] = constel.get_point(point);
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
            static int sps1 = 2;
            static int sps2 = 7;
            static int delay1 = 16;
            static int delay2 = 16;

            // double* taps = new double[sps1*delay1*2+1];
            // size_t len = rrcos(sps1, delay1, 0.35, taps);
            // double* taps2 = new double[2*sps2*delay2+1];
            // size_t len2 = rrcos(sps2, delay2, 0.35, taps2);
            auto prot = root_nyquist(sps1, sps1, 1.0, 0.35, 2*sps1*delay1+1);
            auto prot2 = root_nyquist(sps2, sps2, 1.0, 0.35, 2*sps2*delay2+1);
            fir_resampler<double, std::complex<double>> resamp(sps1, 1, prot.data(), prot.size());
            fir_resampler<double, std::complex<double>> resamp2(1, sps2, prot2.data(), prot2.size());
            // delete[] taps;
            // delete[] taps2;

            size_t n;
            auto h = resamp.operate_n(outs, points, n);
            double* fft = new double[n];
            double* freq = new double[n];
            compute_psd(h, n, fft, freq);

            size_t k;
            auto j = resamp2.operate_n(h, n, k);
            double* real = new double[k];
            double* imag = new double[k];
            for(size_t i = 0; i < k; ++i){
                real[i] = j[i].real();
                imag[i] = j[i].imag();
            }
        
            ImGui::SliderInt("Interp SPS", &sps1, 1, 32);
            ImGui::SliderInt("Decim SPS", &sps2, 1, 32);
            ImGui::SliderInt("Interp Delay", &delay1, 1, 64);
            ImGui::SliderInt("Decim Delay", &delay2, 1, 64);

            if(ImPlot::BeginSubplots("Data", 1, 2, ImVec2(-1, 750))){
            if(ImPlot::BeginPlot("Spectrum")){
                ImPlot::PlotLine("", freq, fft, n);

                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Time Data")){
                ImPlot::PlotScatter("", real, imag, k);
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

