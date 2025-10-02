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
#define RFFT_IMPLEMENTATION

int main(){
    noise<double> sigGen{};
    constellation_qpsk constel{};
    std::vector<double> tx_f = root_nyquist(32, 32, 1.0, 0.35, 11*8*32);
    polyphase_upsampler arb(8, tx_f, 32);

    std::vector<std::complex<double>> iq_data(1024);
    for(size_t i = 0; i < 1024; ++i) iq_data[i] = constel.get_point(sigGen.randomValue(constel.get_bps()));

    auto tx_output = arb.filterN(iq_data);

    fft_transform_radix2(tx_output, false);
    auto h = make_psd(tx_output);

    if(DO_WINDOW){
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
        ImGui::BeginTabBar("Main Tabs");
        if(ImGui::BeginTabItem("FFT")){
            if(ImPlot::BeginPlot("FFT Plot 1", ImVec2(-1, 800))){
                ImPlot::PlotLine("Test Data Filtered", h.data(), h.size());
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

