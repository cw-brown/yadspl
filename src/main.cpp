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
// #include "filter.hpp"

#include "fft.hpp"

void key_call(GLFWwindow* window, int key, int, int action, int){
    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
}

#define DO_WINDOW true
#define RFFT_IMPLEMENTATION

int main(){
    constellation_qpsk QPSK{};
    noise<double> sigGen{};

    std::vector<std::complex<double>> data(1024);
    std::vector<double> realpart(1024);
    std::vector<double> imagpart(1024);
    for(size_t i = 0; i < 1024; ++i){
        data[i] = QPSK.get_point(sigGen.randomValue(QPSK.get_bps()));
        realpart[i] = data[i].real();
        imagpart[i] = data[i].imag();
    }

    auto v = root_nyquist(32, 32, 1.0, 0.35, 11*32*8);
    polyphase_upsampler arb(8, v, 32);

    auto output = arb.filterN(data.data(), 1024);

    fft_transform_radix2(output, false);
    fft_transform_radix2(data, false);

    std::vector<double> plotData = make_psd(output);
    auto minData = make_psd(data);


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
            if(ImPlot::BeginPlot("FFT Plot 1", ImVec2(-1, 750))){
                ImPlot::PlotLine("Test Data Filtered", plotData.data(), plotData.size());
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("FFT Plot 2", ImVec2(-1, 750))){
                ImPlot::PlotLine("Test Data", minData.data(), minData.size());
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("IQ Plots")){
            if(ImPlot::BeginPlot("IQ", ImVec2(-1, 750))){
                ImPlot::PlotScatter("QPSK", realpart.data(), imagpart.data(), 1024);
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

