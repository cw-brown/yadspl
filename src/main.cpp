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
#include "modem.hpp"

#include "fft.hpp"

#define DO_WINDOW true

int main(){
    size_t sps = 4;
    size_t n_filt = 32;

    auto prototype = root_nyquist(n_filt, n_filt, 1.0, 0.35, sps*n_filt*32);

    resampler<double> arb(sps, prototype.data(), prototype.size(), n_filt);
    // auto bank = arb.get_bank();
    // auto deriv = arb.get_deriv_bank();

    double fc = 16;
    double fs = 1024;
    double del_time = 8;
    int num = fs/del_time;
    double* stream = new double[num];
    double* t = new double[num];
    for(int i = 0; i < num; ++i){
        t[i] = i / fs;
        stream[i] = std::sin(2*3.1415926535*fc*t[i]);
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
        if(ImGui::BeginTabItem("Items")){
            if(ImPlot::BeginPlot("Plot 1", ImVec2(-1, 750))){
                ImPlot::PlotStems("", t, stream, num);
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Plot 2", ImVec2(-1, 750))){
                // ImPlot::PlotStems("", upsampled.data(), upsampled.size());
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

