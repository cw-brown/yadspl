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
#include "polynomial.hpp"
#include "noise.hpp"
#include "yadpsl_math.hpp"

void key_call(GLFWwindow* window, int key, int, int action, int){
    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
}

#define DO_WINDOW true

int main(){


    if(DO_WINDOW){
    GLFWwindow* window = glfw_makeNewWindow(1920, 1080, "Yet Another DSP Library", true, true, true);
    ImPlot::CreateContext();
    glfwSetKeyCallback(window, key_call);

    size_t eval_points = 5000;

    static noise<double> gauss;
    // double* signal = gauss.randomSignal(eval_points);
    double* signal = new double[eval_points];
    double* x = new double[eval_points];
    for(size_t i = 0; i < eval_points; ++i){
        x[i] = i;
        if(i >= 2000 && i <= 3000){
            signal[i] = 1.0;
        } else signal[i] = 0.0;
    }


    while(!glfwWindowShouldClose(window)){
        glfwPollEvents();
        glfw_frame();

        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x, ImGui::GetIO().DisplaySize.y));
        ImGuiWindowFlags topbarflags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoCollapse;

        ImGui::Begin("Plottings", (bool*)true, topbarflags);
        ImGui::BeginTabBar("Main Tabs");
        if(ImGui::BeginTabItem("Testing")){
            static fir_filter<double, fir_type::Lowpass_Ideal, window_type::Hamming> windowed_lowpass;
            static int order = 21;
            static float fc = 1.0;
            windowed_lowpass.design(order, fc);
            auto resp = windowed_lowpass.response(eval_points);
            ImGui::SliderInt("Order", &order, 2, 251);
            ImGui::SliderFloat("Cutoff Frequency", &fc, 0.0, 3.14159);
            if(ImPlot::BeginPlot("Impulse Response")){
                ImPlot::SetupAxes("Sample Number", "Amplitude", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::PlotBars("Windowed Lowpass", windowed_lowpass.getTaps(), windowed_lowpass.getNumTaps());
                ImPlot::EndPlot();
            } if(ImPlot::BeginPlot("Frequency Response")){
                ImPlot::SetupAxes("Frequency (rad/s)", "Magnitude (dB)", ImPlotAxisFlags_AutoFit, 0);
                ImPlot::SetupAxisLimits(ImAxis_Y1, -80, 10);
                ImPlot::PlotLine("Windowed Lowpass", resp.first, resp.second, eval_points);
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("Noise Sample")){
            if(ImPlot::BeginPlot("Noise")){
                ImPlot::PlotLine("Sample", x, signal, eval_points);
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("Filtered Sample")){
            static fir_filter<double, fir_type::Lowpass_Ideal> lowpass;
            static int order = 51;
            static float fc = 1.5;
            lowpass.design(order, fc);
            lowpass.synthesize();
            ImGui::SliderInt("Filter Order", &order, 3, 251);
            ImGui::SliderFloat("Cutoff Frequency", &fc, 0.0, 3.141592);
            if(ImPlot::BeginPlot("Filtered")){
                ImPlot::PlotLine("Sample", x, lowpass.filter(signal, eval_points), eval_points);
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("DFT Of Signal")){

            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("Bandpass Filter")){

            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("Bandstop Filter")){
            
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

