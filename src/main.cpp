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
// #include "constellations2.hpp"

// #include "filter.hpp"

void key_call(GLFWwindow* window, int key, int, int action, int){
    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
}

#define DO_WINDOW true

int main(){
    constellation_16qam constel{};

    symbol_recovery rec(8, 3.141/50.0, 32, 1.5, constel, 0.35);
    auto bank = rec.get_bank();

    noise<double> sigGen;

    std::ring<double> realPart(1024);
    std::ring<double> imagPart(1024);

    if(DO_WINDOW){
    GLFWwindow* window = glfw_makeNewWindow(1920, 1080, "Yet Another DSP Library", true, true, true);
    ImPlot::CreateContext();
    glfwSetKeyCallback(window, key_call);


    while(!glfwWindowShouldClose(window)){
        glfw_frame();

        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x, ImGui::GetIO().DisplaySize.y));
        ImGuiWindowFlags topbarflags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoCollapse;
        ImPlotAxisFlags axisflags = ImPlotAxisFlags_AutoFit;

        // Start to make data
        auto point = sigGen.randomConstellationPoint(constel);
        realPart.push_back(point.real());
        imagPart.push_back(point.imag());
        

        ImGui::Begin("Plottings", nullptr, topbarflags);
        ImGui::BeginTabBar("Main Tabs");
        if(ImGui::BeginTabItem("Test Data")){
            if(ImPlot::BeginPlot("Testing Data", ImVec2(-1, 700))){
                ImPlot::SetupAxesLimits(-1.5, 1.5, -1.5, 1.5);
                ImPlot::PlotScatter("Test IQ", realPart.data(), imagPart.data(), realPart.size());
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("Test Feature")){
            if(ImPlot::BeginPlot("Stem", ImVec2(-1, 700))){
                ImPlot::SetupAxis(ImAxis_X1, "Sample", axisflags);
                ImPlot::SetupAxis(ImAxis_Y1, "Amplitude", axisflags);
                ImPlot::PlotBars("FIR Filter", bank.get_prototype().data(), bank.get_prototype().size());
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("Split Filter")){
            static int arm = 0;
            ImGui::SliderInt("Arm Number", &arm, 0, bank.get_size()-1);
            auto s = bank.get_arm(arm);
            if(ImPlot::BeginPlot("Normal", ImVec2(-1, 700))){
                ImPlot::SetupAxis(ImAxis_X1, "Sample", axisflags);
                ImPlot::SetupAxis(ImAxis_Y1, "Amplitude", axisflags);
                ImPlot::PlotBars("Arm", s.taps().data(), s.size());
                ImPlot::EndPlot();
            }
            auto d = bank.get_deriv_arm(arm);
            if(ImPlot::BeginPlot("Derivative", ImVec2(-1, 700))){
                ImPlot::SetupAxis(ImAxis_X1, "Sample", axisflags);
                ImPlot::SetupAxis(ImAxis_Y1, "Amplitude", axisflags);
                ImPlot::PlotBars("Arm", d.taps().data(), d.size());
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

